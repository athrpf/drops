/// \file levelset.cpp
/// \brief levelset equation for two phase flow problems
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include "levelset/levelset.h"
#include "levelset/fastmarch.h"
#include "num/lattice-eval.h"
#include "num/quadrature.h"
#include <fstream>

namespace DROPS
{

inline double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha);
}

/// \brief Base class for all surface tension accumulators
class SurfTensAccumulatorCL : public TetraAccumulatorCL
{
  protected:
    VecDescCL  SmPhi_;
    const BndDataCL<>& lsetbnd_;
    VecDescCL& f;
    SMatrixCL<3,3> T;
    InterfaceTriangleCL triangle;

  public:
    SurfTensAccumulatorCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma)
     : SmPhi_(ls.Phi), lsetbnd_(ls.GetBndData()), f(f_Gamma)
    { ls.MaybeSmooth( SmPhi_.Data); }

    void begin_accumulation ()
    {
        // uncomment for Geomview output
        //std::ofstream fil("surf.off");
        //fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";
    }
    void finalize_accumulation()
    {
        // uncomment for Geomview output
        //fil << "}\n";
    }
    ///\brief Do setup of f_Gamma on given tetra
    virtual void visit (const TetraCL&)= 0;

};

/// \brief Accumulator for the (artificial) constant surface force, mainly used for numerical test cases.
///
/// Computes the integral
///         \f[ \sigma \int_\Gamma v \textbf n ds. \f]
class ConstSurfTensAccumulatorCL : public SurfTensAccumulatorCL
{
  private:
    const double sigma_;
    const RefRuleCL RegRef_;

    IdxT Numb[14];
    LocalP2CL<> p1abs_p[4][8], p1abs_n[4][8]; // extended basis functions on pos./neg. part, resp., for each of the 8 regular children
    LocalP2CL<> loc_phi;

  public:
    ConstSurfTensAccumulatorCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sigma_(sigma), RegRef_(GetRefRule( RegRefRuleC)) {}

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new ConstSurfTensAccumulatorCL ( *this); };
};

void ConstSurfTensAccumulatorCL::visit (const TetraCL& t)
{
    const Uint idx_f= f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();

    loc_phi.assign( t, SmPhi_, lsetbnd_);
    triangle.Init( t, loc_phi);
    if (!triangle.Intersects())
        return;

    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
    }
    for (int xv=0; xv<4; ++xv)
    { // collect data on all extended DoF
        Numb[xv+10]= velXfem && Numb[xv]!=NoIdx ? f.RowIdx->GetXidx()[Numb[xv]] : NoIdx;
    }

    if (velXfem)
        P2RidgeDiscCL::GetExtBasisOnChildren( p1abs_p, p1abs_n, loc_phi);

    for (int ch=0; ch<8; ++ch)
    {
        if (!triangle.ComputeForChild(ch)) // no patch for this child
            continue;

//patch.WriteGeom( fil);

        BaryCoordCL BaryPQR, BarySQR;
        for (int i=0; i<3; ++i)
        {
            // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
            BaryPQR+= triangle.GetBary(i);
            BarySQR+= triangle.GetBary(i+1);
        }
        BaryPQR/= 3; BarySQR/= 3;

        Point3DCL n; // senkrecht auf PQ, PR, nach aussen gerichtet...
        cross_product( n, triangle.GetPoint(1)-triangle.GetPoint(0), triangle.GetPoint(2)-triangle.GetPoint(0));
        n/= n.norm();

        const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
        const int find_sign= triangle.GetNumSign( 1) ? 1 : -1;
        Point3DCL pos_dir;
        for (Uint i=0; i<4; ++i) // compute vector with positive direction
        {
            const Uint vert= data.Vertices[i];
            if (triangle.GetSign(vert)==find_sign)
            {
                const Point3DCL signedPoint= vert<4 ? t.GetVertex(vert)->GetCoord() : GetBaryCenter( *t.GetEdge(vert-4));
                pos_dir= signedPoint - triangle.GetPoint(0);
                if (find_sign == -1) pos_dir= -pos_dir;
                break;
            }
        }
        if (inner_prod( n, pos_dir) < 0) n= -n;

        double val_hat[4];
        for (int v=0; v<14; ++v)
        {
            if (Numb[v]==NoIdx) continue;

            for (Uint k=0; k<triangle.GetNumPoints(); ++k)
                // values of basis function in P,Q,R,S. Note: p1abs_p==p1abs_n on \f$Gamma_h\f$
                val_hat[k]= v<10 ? FE_P2CL::H(v,triangle.GetBary(k)) : p1abs_p[v-10][ch](triangle.GetBary(k));

            double v_Bary= v<10 ? FE_P2CL::H(v,BaryPQR) : p1abs_p[v-10][ch](BaryPQR),
                sum_v= 0;
            for (int k=0; k<3; ++k)
                 sum_v+= val_hat[k];

            if (triangle.IsQuadrilateral())
            {
                double sum_vSQR= 0;
                for (int k=1; k<4; ++k)
                    sum_vSQR+= val_hat[k];
                sum_v+= triangle.GetAreaFrac() * sum_vSQR;
                v_Bary+= triangle.GetAreaFrac() * (v<10 ? FE_P2CL::H(v,BarySQR) : p1abs_p[v-10][ch](BarySQR));
            }

            // Quadraturformel auf Dreieck, exakt bis zum Grad 2
            const double int_v= (1./12)*sum_v + 0.75 * v_Bary;
            const double C= sigma_*triangle.GetAbsDet()/2;
            for (int i=0; i<3; ++i)
                f.Data[Numb[v]+i]-= C * int_v*n[i];
        }
    } // Ende der for-Schleife ueber die Kinder
}


/// \brief Accumulator for the naive Laplace-Beltrami discretization of the CSF term.
///
/// Computes the integral
///         \f[ \sigma \int_\Gamma \kappa v \textbf n_h ds = \sigma \int_\Gamma P_h\nabla id P_h\nabla v ds \f]
/// with \f$P_h = I - \textbf n_h\cdot \textbf n_h^T\f$.
/// Discretization error in \f$ H^1(\Omega)\f$ dual norm has only order 1/2 w.r.t. the grid size at the interface.
/// Better use the improved Laplace-Beltrami discretization with better discretization order.
class NaiveLaplaceBeltramiAccuCL : public SurfTensAccumulatorCL
{
  private:
    const double sigma_;

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    IdxT Numb[10];

  public:
    NaiveLaplaceBeltramiAccuCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sigma_(sigma)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new NaiveLaplaceBeltramiAccuCL ( *this); };
};

void NaiveLaplaceBeltramiAccuCL::visit( const TetraCL& t)
// computes the integral
//         sigma \int_\Gamma \kappa v n ds = sigma \int_\Gamma grad id grad v ds
{
    const Uint  idx_f= f.RowIdx->GetIdx();
    double det;

    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
    }

    triangle.Init( t, SmPhi_, lsetbnd_);

    for (int ch=0; ch<8; ++ch)
    {
        if (!triangle.ComputeForChild(ch)) // no patch for this child
            continue;

//patch.WriteGeom( fil);

        BaryCoordCL BaryPQR, BarySQR;
        for (int i=0; i<3; ++i)
        {
            // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
            BaryPQR+= triangle.GetBary(i);
            BarySQR+= triangle.GetBary(i+1);
        }

        const double C= triangle.GetAbsDet()/6*sigma_;  // 1/6 for quad. rule

        for (int v=0; v<10; ++v)
        {
            if (Numb[v]==NoIdx) continue;

            LocalP1CL<Point3DCL> gradv; // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
            for (int node=0; node<4; ++node)
                gradv[node]= Grad[v][node];

            Point3DCL gr= gradv( BaryPQR); // gr= grad v(P) + grad v(Q) + grad v(R)

            if (triangle.IsQuadrilateral())
                gr+= triangle.GetAreaFrac() * gradv( BarySQR);
            // nun gilt:
            // gr = [grad v(P)+...] + (a+b-1)[grad v(S)+...]

            for (int i=0; i<3; ++i)
            {
                const double val= inner_prod( gr, triangle.GetGradId(i));
                f.Data[Numb[v]+i]-= C *val;
            }
        }
    } // Ende der for-Schleife ueber die Kinder
}

/// \brief Accumulator for the improved Laplace-Beltrami discretization of the CSF term.
///
/// Computes the integral
///         \f[ \sigma \int_\Gamma \kappa v \textbf n_h ds = \sigma \int_\Gamma \hat P_h \nabla id \hat P_h\nabla v ds \f]
/// with \f$\hat P_h = \tilde P_h P_h, P_h = I - \textbf n_h\cdot \textbf n_h^T, \tilde P_h = I - \tilde\textbf n_h\cdot \tilde\textbf n_h^T \f$.
/// Discretization error in \f$ H^1(\Omega)\f$ dual norm has order 1 w.r.t. the grid size at the interface, numerical experiments even indicate order 1.5.
class ImprovedLaplaceBeltramiAccuCL : public SurfTensAccumulatorCL
{
  private:
    const double sigma_;

    LocalP1CL<Point3DCL> Grad[10], GradRef[10];
    IdxT Numb[10];
    LocalP2CL<> velR_p[4][8], velR_n[4][8]; // for P2R basis on children
    LocalP2CL<> loc_phi;

  public:
    ImprovedLaplaceBeltramiAccuCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, double sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sigma_(sigma)
    { P2DiscCL::GetGradientsOnRef( GradRef); }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new ImprovedLaplaceBeltramiAccuCL ( *this); };
};

void ImprovedLaplaceBeltramiAccuCL::visit ( const TetraCL& t)
{
    const Uint idx_f=   f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();
    double det;

    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
    LocalP1CL<Point3DCL> n;

    loc_phi.assign( t, SmPhi_, lsetbnd_);
    triangle.Init( t, loc_phi);
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
        for (int k=0; k<4; ++k)
            n[k]+= triangle.GetPhi(v)*Grad[v][k];
    }

    for (int ch=0; ch<8; ++ch)
    {
        if (!triangle.ComputeForChild(ch)) // no patch for this child
            continue;

//patch.WriteGeom( fil);
        BaryCoordCL BaryPQR, BarySQR;
        for (int i=0; i<3; ++i)
        {
            // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
            BaryPQR+= triangle.GetBary(i);
            BarySQR+= triangle.GetBary(i+1);
        }
        BaryPQR/= 3.;    BarySQR/= 3.;

        typedef SArrayCL<Point3DCL,3> ProjT;
        GridFunctionCL<ProjT> GradId( ProjT(), 6);  // values in P, Q, R, S, BaryPQR, BarySQR
        for (int p=0; p<6; ++p)
        {
            Point3DCL np= n( p<4 ? triangle.GetBary(p) : p==4 ? BaryPQR : BarySQR);
            if (np.norm()>1e-8) np/= np.norm();
            for (int i=0; i<3; ++i)
                GradId[p][i]= triangle.ApplyProj( std_basis<3>(i+1) - np[i]*np);
            //                     GradId[p][i]= std_basis<3>(i+1) - np[i]*np;
        }
        const double C= triangle.GetAbsDet()*sigma_/2.;
        if (velXfem)
            P2RidgeDiscCL::GetExtBasisOnChildren( velR_p, velR_n, loc_phi);
        for (int v=0; v<(velXfem ? 14 : 10); ++v)
        {
            const IdxT Numbv= v<10 ? Numb[v] : (velXfem && Numb[v-10]!=NoIdx ? f.RowIdx->GetXidx()[Numb[v-10]] : NoIdx);
            if (Numbv==NoIdx) continue;

            LocalP1CL<Point3DCL> gradv; // gradv = gradient of hat function for dof v
            if (v<10) // std basis function
                for (int node=0; node<4; ++node)
                    gradv[node]= Grad[v][node];
            else // extended basis function: tangential derivative is the same for pos./neg. part, ie., P_h grad(vx_p) == P_h grad(vx_n). W.l.o.g. take pos. part for computation.
                P2DiscCL::GetFuncGradient( gradv, velR_p[v-10][ch], Grad);

            for (int i=0; i<3; ++i)
            {
                double intSum= 0; // sum of the integrand in PQR, SQR
                for (int k=0; k<3; ++k)
                {
                    intSum+= inner_prod( GradId[k][i], gradv(triangle.GetBary(k)));
                    if (triangle.IsQuadrilateral())
                        intSum+= triangle.GetAreaFrac() * inner_prod( GradId[k+1][i], gradv(triangle.GetBary(k+1)));
                }
                double intBary= inner_prod( GradId[4][i], gradv(BaryPQR));
                if (triangle.IsQuadrilateral())
                    intBary+= triangle.GetAreaFrac() * inner_prod( GradId[5][i], gradv(BarySQR));
                f.Data[Numbv+i]-= C *(intSum/12. + 0.75*intBary);
            }
        }
    } // Ende der for-Schleife ueber die Kinder
}

void SF_ImprovedLaplBeltramiOnTriangle( const TetraCL& t, const BaryCoordCL * const p,
                                        const InterfaceTriangleCL&  triangle, const LocalP1CL<Point3DCL> Grad_f[10], const IdxT Numb[10],
                                        instat_scalar_fun_ptr sigma, const Quad5_2DCL<Point3DCL> e[3],
                                        double det, VectorCL& f)
{
    static Quad5_2DCL<Point3DCL> Grad[10]; // Gradients of the P2-basis-functions
    Quad5_2DCL<Point3DCL> n;

    for (int v=0; v<10; ++v)
    {
        Grad[v].assign( Grad_f[v], p);
        n+= triangle.GetPhi(v)*Grad[v];
    }
    for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();

    Quad5_2DCL<> qsigma( t, p, sigma),  // surface tension
                 q1;                    // Term 1

    Quad5_2DCL<Point3DCL> qPhPhte,      // Common term in Term 1 and Term 2
                          qsigmaPhPhte; // for Term 1

    for (int i= 0; i < 3; ++i)
    {
        qPhPhte= (e[i] - dot(e[i],n)*n);
        qPhPhte.apply( triangle, &InterfaceTriangleCL::ApplyProj);
        qsigmaPhPhte= qsigma*qPhPhte;
        for (int v= 0; v < 10; ++v)
        {
            if (Numb[v]==NoIdx) continue;
            q1= dot (qsigmaPhPhte, Grad[v]);
            f[Numb[v]+i]-= q1.quad( det);
        }
    }
}

void SF_ImprovedLaplBeltramiOnTriangle( const TetraCL& t, const BaryCoordCL * const p,
    const InterfaceTriangleCL&  triangle, const LocalP1CL<Point3DCL> Grad_f[10], const IdxT Numb[10],
    const Quad5_2DCL<Point3DCL> e[3], double det, VectorCL& f, const SurfaceTensionCL& sf)
{
    if (sf.GetInputMethod() == Sigma_X)
    {
        SF_ImprovedLaplBeltramiOnTriangle( t, p, triangle, Grad_f, Numb, sf.GetSigma(), e, det, f);
        return;
    }
    static Quad5_2DCL<>          p2[10];   // P2-Hat-Functions...
    static Quad5_2DCL<Point3DCL> Grad[10]; // and their gradients
    Quad5_2DCL<Point3DCL> n;
    P2DiscCL::GetP2Basis( p2, p);
    for (int v=0; v<10; ++v)
    {
        Grad[v].assign( Grad_f[v], p);
        n+= triangle.GetPhi(v)*Grad[v];
    }
    for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();

    Quad5_2DCL<> qsigma, q1;
    Quad5_2DCL<Point3DCL> qPhPhte,                         // Common term in Term 1 and Term 2
                          qsigmaPhPhte;                    // for Term 1
    sf.ComputeSF(t, p, qsigma);

    for (int i= 0; i < 3; ++i)
    {
        qPhPhte= (e[i] - dot(e[i],n)*n);
        qPhPhte.apply( triangle, &InterfaceTriangleCL::ApplyProj);
        qsigmaPhPhte= qsigma*qPhPhte;
        for (int v= 0; v < 10; ++v)
        {
            if (Numb[v]==NoIdx) continue;
            q1= dot (qsigmaPhPhte, Grad[v]);
            f[Numb[v]+i]-= q1.quad( det);
        }
    }
}

/// \brief Accumulator for the improved Laplace-Beltrami discretization of the CSF term with variable surface tension coefficient.
class VarImprovedLaplaceBeltramiAccuCL : public SurfTensAccumulatorCL
{
  private:
    const SurfaceTensionCL& sf_;

    LocalP1CL<Point3DCL> GradRef[10], Grad[10];
    IdxT Numb[10];
    Quad5_2DCL<Point3DCL> e[3];

  public:
    VarImprovedLaplaceBeltramiAccuCL( const LevelsetP2CL& ls, VecDescCL& f_Gamma, const SurfaceTensionCL& sigma)
     : SurfTensAccumulatorCL( ls, f_Gamma), sf_(sigma)
    {
        P2DiscCL::GetGradientsOnRef( GradRef);
        for (int i= 0; i<3; ++i)
            e[i]= std_basis<3>( i + 1);
    }

    void visit (const TetraCL&);

    TetraAccumulatorCL* clone (int /*tid*/) { return new VarImprovedLaplaceBeltramiAccuCL ( *this); };
};

void VarImprovedLaplaceBeltramiAccuCL::visit( const TetraCL& t)
// computes the integral sigma \int_\Gamma \kappa v n ds = sigma
// \int_\Gamma grad id grad v ds
{
    const Uint idx_f= f.RowIdx->GetIdx();
    double det;

    triangle.Init( t, SmPhi_, lsetbnd_);

    for (int v= 0; v < 10; ++v)
    { // collect data on all DoF
        const UnknownHandleCL& unk= v<4 ? t.GetVertex(v)->Unknowns : t.GetEdge(v-4)->Unknowns;
        Numb[v]= unk.Exist( idx_f) ? unk( idx_f) : NoIdx;
    }
    GetTrafoTr( T, det, t);
    P2DiscCL::GetGradients( Grad, GradRef, T);

    for (int ch= 0; ch < 8; ++ch)
    {
        triangle.ComputeForChild( ch);
        for (int tri= 0; tri < triangle.GetNumTriangles(); ++tri)
            SF_ImprovedLaplBeltramiOnTriangle( t, &triangle.GetBary( tri),
                    triangle, Grad,  Numb, e, triangle.GetAbsDet( tri), f.Data, sf_);
    } // Ende der for-Schleife ueber die Kinder
}

void MarkInterface ( scalar_fun_ptr DistFct, double width, MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, it)
    {
        double d= 1e99;
        int num_pos= 0;
        for (Uint j=0; j<10; ++j)
        {
            const double dist= j<4 ? DistFct( it->GetVertex( j)->GetCoord())
                                   : DistFct( GetBaryCenter( *it->GetEdge(j-4)));
            if (dist>=0) ++num_pos;
            d= std::min( d, std::abs( dist));
        }

        const bool vzw= num_pos!=0 && num_pos!=10; // change of sign
        if (d<=width || vzw)
            it->SetRegRefMark();
    }
}

void MarkInterface ( const LevelsetP2CL::const_DiscSolCL& lset, double width, MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/-1, it)
    {
        double d= 1e99;
        int num_pos= 0;
        for (Uint j=0; j<10; ++j)
        {
            const double dist= j<4 ? lset.val( *it->GetVertex( j))
                                   : lset.val( *it->GetEdge(j-4));
            if (dist>=0) ++num_pos;
            d= std::min( d, std::abs( dist));
        }

        const bool vzw= num_pos!=0 && num_pos!=10; // change of sign
        if (d<=width || vzw)
            it->SetRegRefMark();
    }
}


//*****************************************************************************
//                               LevelsetP2CL
//*****************************************************************************

void LevelsetP2CL::Init( scalar_fun_ptr phi0)
{
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    for (MultiGridCL::TriangVertexIteratorCL it= MG_.GetTriangVertexBegin(lvl),
        end= MG_.GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        if ( it->Unknowns.Exist(idx))
        Phi.Data[it->Unknowns(idx)]= phi0( it->GetCoord());
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(lvl),
        end= MG_.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        if ( it->Unknowns.Exist(idx))
        Phi.Data[it->Unknowns(idx)]= phi0( GetBaryCenter( *it));
    }
}


void LevelsetP2CL::CreateNumbering( Uint level, IdxDescCL* idx, match_fun match)
{
    idx->CreateNumbering( level, MG_, BndData_, match);
}


void LevelsetP2CL::Reparam( int method, bool Periodic)
/** \param method How to perform the reparametrization (see description of ReparamFactoryCL for details)
    \param Periodic: If true, a special variant of the algorithm for periodic boundaries is used.
*/
{
    std::auto_ptr<ReparamCL> reparam= ReparamFactoryCL::GetReparam( MG_, Phi, method, Periodic, &BndData_, perDirections);
    reparam->Perform();
}

void LevelsetP2CL::AccumulateBndIntegral( VecDescCL& f) const
{
    SurfTensAccumulatorCL* accu;

    switch (SF_)
    {
      case SF_LB:
          accu= new NaiveLaplaceBeltramiAccuCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.)); break;
      case SF_Const:
          accu= new ConstSurfTensAccumulatorCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.)); break;
      case SF_ImprovedLB:
          accu= new ImprovedLaplaceBeltramiAccuCL( *this, f, sf_.GetSigma()(std_basis<3>(0), 0.)); break;
      case SF_ImprovedLBVar:
          accu= new VarImprovedLaplaceBeltramiAccuCL( *this, f, sf_); break;
      default:
        throw DROPSErrCL("LevelsetP2CL::AccumulateBndIntegral not implemented for this SurfaceForceT");
    }
    TetraAccumulatorTupleCL accus;
    accus.push_back( accu);
    accumulate( accus, MG_, Phi.RowIdx->TriangLevel(), Phi.RowIdx->GetMatchingFunction(), Phi.RowIdx->GetBndInfo());

    delete accu;
}

double LevelsetP2CL::GetVolume( double translation, int l) const
{
    if (l==0)
        ++l;
    double Volume= l > 0 ? GetVolume_Composite( translation, l)
                         : GetVolume_Extrapolation( translation, -l);

#ifdef _PAR
    Volume = ProcCL::GlobalSum(Volume);
#endif
    return Volume;
}

double LevelsetP2CL::GetVolume_Extrapolation( double translation, int l) const
{
    double vol = 0.;
    QuadDomainCL qdom;
    DROPS::ExtrapolationToZeroCL extra( l, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( l, DROPS::HarmonicSubdivisionCL());
    LocalP2CL<> loc_phi;
    DROPS_FOR_TRIANG_TETRA( MG_, idx.TriangLevel(), it) {
        loc_phi.assign(*it,Phi,GetBndData());
        loc_phi+= translation;
        make_ExtrapolatedQuad5Domain( qdom, loc_phi, extra);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        vol+=quad( integrand, it->GetVolume()*6., qdom, NegTetraC);
    }
    return vol;
}

double LevelsetP2CL::GetVolume_Composite( double translation, int l) const
{
    const PrincipalLatticeCL& lat= PrincipalLatticeCL::instance ( l);
    double vol = 0.;
    std::valarray<double> ls_values (lat.vertex_size());
    QuadDomainCL qdom;
    LocalP2CL<> loc_phi;
    TetraPartitionCL partition;
    DROPS_FOR_TRIANG_TETRA( MG_, idx.TriangLevel(), it) {
        loc_phi.assign(*it,Phi,GetBndData());
        loc_phi+= translation;
        evaluate_on_vertexes (loc_phi, lat, Addr(ls_values));
        partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL>(lat, ls_values);
        make_CompositeQuad5Domain( qdom, partition);
        DROPS::GridFunctionCL<> integrand( 1., qdom.vertex_size());
        vol+=quad( integrand, it->GetVolume()*6., qdom, NegTetraC);
    }
    return vol;
}

double LevelsetP2CL::AdjustVolume (double vol, double tol, double surface, int l) const
{
    tol*=vol;

    double v0=GetVolume(0., l)-vol;
    if (std::abs(v0)<=tol) return 0;

    double d0=0, d1=v0*(surface != 0. ? 1.1/surface : 0.23/std::pow(vol,2./3.));
    // Hinweis: surf(Kugel) = [3/4/pi*vol(Kugel)]^(2/3) * 4pi
    double v1=GetVolume(d1, l)-vol;
    if (std::abs(v1)<=tol) return d1;

    // Sekantenverfahren fuer Startwert
    while (v1*v0 > 0) // gleiches Vorzeichen
    {
        const double d2=d1-1.2*v1*(d1-d0)/(v1-v0);
        d0=d1; d1=d2; v0=v1; v1=GetVolume(d1, l)-vol;
        if (std::abs(v1)<=tol) return d1;
    }

    // Anderson-Bjoerk fuer genauen Wert
    while (true)
    {
        const double d2=(v1*d0-v0*d1)/(v1-v0),
                     v2=GetVolume(d2,l)-vol;
        if (std::abs(v2)<=tol) return d2;

        if (v2*v1 < 0) // ungleiches Vorzeichen
          { d0=d1; d1=d2; v0=v1; v1=v2; }
        else
          { const double c=1.0-v2/v1; d1=d2; v1=v2; v0*= c>0 ? c : 0.5; }
    }
}

void LevelsetP2CL::SmoothPhi( VectorCL& SmPhi, double diff) const
{
    Comment("Smoothing for curvature calculation\n", DebugDiscretizeC);
    MatrixCL M, A, C;
    SetupSmoothSystem( M, A);
    C.LinComb( 1, M, diff, A);
#ifndef _PAR
    SSORPcCL pc;
    PCG_SsorCL pcg( pc, 500, 1e-10);
    pcg.Solve( C, SmPhi, M*Phi.Data);
    __UNUSED__ double inf_norm= supnorm( SmPhi-Phi.Data);
#else
    ParJac0CL  JACPc (idx);
    typedef ParPCGSolverCL<ParJac0CL> JacPCGSolverT;
    JacPCGSolverT cg( 500, 1e-10, idx, JACPc);
    cg.Solve( C, SmPhi, M*Phi.Data);
    __UNUSED__ const double inf_norm= ProcCL::GlobalMax(supnorm( SmPhi-Phi.Data));
#endif
    Comment("||SmPhi - Phi||_oo = " <<inf_norm<< std::endl, DebugDiscretizeC);
}

void LevelsetP2CL::SetupSmoothSystem( MatrixCL& M, MatrixCL& A) const
// used for smoothing of Phi before computing curvature term
//
// M = mass matrix for P2 elements
// A = stiffness matrix for P2 elements
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns();
    const Uint lvl= Phi.GetLevel();

#ifndef _PAR
    __UNUSED__ const IdxT allnum_unks= num_unks;
#else
    __UNUSED__ const IdxT allnum_unks= ProcCL::GlobalSum(num_unks);
#endif
    Comment("entering Levelset::SetupSmoothSystem: " << allnum_unks << " levelset unknowns.\n", DebugDiscretizeC);

    SparseMatBuilderCL<double> Mb(&M, num_unks, num_unks);
    SparseMatBuilderCL<double> Ab(&A, num_unks, num_unks);

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3>     T;
    P2DiscCL::GetGradientsOnRef( GradRef);

    IdxT         Numb[10];
    double det, absdet;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        GetLocalNumbP2NoBnd( Numb, *sit, *Phi.RowIdx);

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            for(int j=0; j<10; ++j)
            {
                // M_ij = ( v_j, v_i),    A_ij = ( grad v_j, grad v_i)
                Mb( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j)*absdet;
                Ab( Numb[i], Numb[j])+= Quad2CL<>(dot( Grad[j], Grad[i])).quad( absdet);
            }
        }
    }
    Mb.Build();
    Ab.Build();
}

void LevelsetP2CL::GetMaxMinGradPhi(double& maxGradPhi, double& minGradPhi) const
{
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det;
    InterfacePatchCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
    maxGradPhi= -1.;
    minGradPhi= 1e99;

    double maxNorm;
    double minNorm;

    DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), it)
    {
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
        patch.Init( *it, Phi, BndData_);

        // compute maximal norm of grad Phi
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
            gradPhi+= patch.GetPhi(v)*Grad[v];

        VectorCL normGrad( 5);
        for (int v=0; v<5; ++v) // init normGrad
            normGrad[v]= norm( gradPhi[v]);
        maxNorm= normGrad.max();
        minNorm= normGrad.min();
        if (maxNorm > maxGradPhi) maxGradPhi= maxNorm;
        if (minNorm < minGradPhi && patch.Intersects()) minGradPhi= minNorm;
    }
#ifdef _PAR
    maxGradPhi= ProcCL::GlobalMax( maxNorm= maxGradPhi);
    minGradPhi= ProcCL::GlobalMin( minNorm= minGradPhi);
#endif
}

//*****************************************************************************
//                               LevelsetRepairCL
//*****************************************************************************

void LevelsetRepairCL::pre_refine()
{
#ifndef _PAR
    p2repair_= std::auto_ptr<RepairP2CL<double> >(
        new RepairP2CL<double>( ls_.GetMG(), ls_.Phi, ls_.GetBndData()));
#else
    /// Tell parallel multigrid about the location of the DOF
    GetPMG().AttachTo( &ls_.Phi, &ls_.GetBndData());
#endif
}

void
LevelsetRepairCL::post_refine ()
/// Do all things to complete the repairing of the FE level-set function
{
    VecDescCL loc_phi;
    IdxDescCL loc_lidx( P2_FE);
    VecDescCL& phi= ls_.Phi;
    match_fun match= ls_.GetMG().GetBnd().GetMatchFun();

    ls_.CreateNumbering( ls_.GetMG().GetLastLevel(), &loc_lidx, match);
    loc_phi.SetIdx( &loc_lidx);
#ifndef _PAR
    p2repair_->repair( loc_phi);
#else
    GetPMG().HandleNewIdx(&ls_.idx, &loc_phi);
    RepairAfterRefineP2( ls_.GetSolution( phi), loc_phi);
    GetPMG().CompleteRepair( &loc_phi);
#endif

    phi.Clear( phi.t);
    ls_.DeleteNumbering( phi.RowIdx);
    ls_.idx.swap( loc_lidx);
    phi.SetIdx( &ls_.idx);
    phi.Data= loc_phi.Data;
}

} // end of namespace DROPS


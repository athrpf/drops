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
#include <fstream>

namespace DROPS
{

inline double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha);
}

void SF_ConstForce( const MultiGridCL& MG, const VecDescCL& SmPhi, const BndDataCL<>& lsetbnd, double sigma, VecDescCL& f)
// computes the integral
//         sigma \int_\Gamma v n ds
// used by levelset/prJump.cpp for testing FE pressure spaces.
{
    const Uint idx_f= f.RowIdx->GetIdx();
    const bool velXfem= f.RowIdx->IsExtended();
    IdxT Numb[14];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    LocalP2CL<> p1abs_p[4][8], p1abs_n[4][8]; // extended basis functions on pos./neg. part, resp., for each of the 8 regular children
    LocalP2CL<> loc_phi;
    SMatrixCL<3,3> T;
    InterfaceTriangleCL triangle;

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    for (MultiGridCL::const_TriangTetraIteratorCL it=MG.GetTriangTetraBegin(), end=MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        loc_phi.assign( *it, SmPhi, lsetbnd);
        triangle.Init( *it, loc_phi);
        if (!triangle.Intersects())
            continue;

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
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

            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            const int find_sign= triangle.GetNumSign( 1) ? 1 : -1;
            Point3DCL pos_dir;
            for (Uint i=0; i<4; ++i) // compute vector with positive direction
            {
                const Uint vert= data.Vertices[i];
                if (triangle.GetSign(vert)==find_sign)
                {
                    const Point3DCL signedPoint= vert<4 ? it->GetVertex(vert)->GetCoord() : GetBaryCenter( *it->GetEdge(vert-4));
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
                const double C= sigma*triangle.GetAbsDet()/2;
                for (int i=0; i<3; ++i)
                    f.Data[Numb[v]+i]-= C * int_v*n[i];
            }
        } // Ende der for-Schleife ueber die Kinder
    }
//fil << "}\n";
}

void SF_LaplBeltrami( const MultiGridCL& MG, const VecDescCL& SmPhi, const BndDataCL<>& lsetbnd, double sigma, VecDescCL& f)
// computes the integral
//         sigma \int_\Gamma \kappa v n ds = sigma \int_\Gamma grad id grad v ds
{
    const Uint  idx_f=     f.RowIdx->GetIdx();
    IdxT        Numb[10];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det;
    InterfaceTriangleCL triangle;

    P2DiscCL::GetGradientsOnRef( GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL it=MG.GetTriangTetraBegin(), end=MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
            Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
        }

        triangle.Init( *it, SmPhi, lsetbnd);

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

            const double C= triangle.GetAbsDet()/6*sigma;  // 1/6 for quad. rule

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
//fil << "}\n";
}

void SF_ImprovedLaplBeltrami( const MultiGridCL& MG, const VecDescCL& SmPhi, const BndDataCL<>& lsetbnd, double sigma, VecDescCL& f)
// computes the integral
//         sigma \int_\Gamma \kappa v n ds = sigma \int_\Gamma grad id grad v ds
{
  const Uint idx_f=   f.RowIdx->GetIdx();
  const bool velXfem= f.RowIdx->IsExtended();
  IdxT       Numb[10];
//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

  LocalP1CL<Point3DCL> Grad[10], GradRef[10];
  SMatrixCL<3,3> T;
  double det;
  InterfaceTriangleCL triangle;
  LocalP2CL<> velR_p[4][8], velR_n[4][8]; // for P2R basis on children
  LocalP2CL<> loc_phi;

  P2DiscCL::GetGradientsOnRef( GradRef);

  for (MultiGridCL::const_TriangTetraIteratorCL it=MG.GetTriangTetraBegin(), end=MG.GetTriangTetraEnd();
       it!=end; ++it)
  {
    GetTrafoTr( T, det, *it);
    P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
    LocalP1CL<Point3DCL> n;

    loc_phi.assign( *it, SmPhi, lsetbnd);
    triangle.Init( *it, loc_phi);
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
      const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
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
      const double C= triangle.GetAbsDet()*sigma/2.;
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
//fil << "}\n";
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
    const Quad5_2DCL<Point3DCL> e[3], double det, VectorCL& f, SurfaceTensionCL& sf)
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

void SF_ImprovedLaplBeltrami( const MultiGridCL& MG, const VecDescCL& SmPhi, const BndDataCL<>& lsetbnd, VecDescCL& f, SurfaceTensionCL& sf)
// computes the integral sigma \int_\Gamma \kappa v n ds = sigma
// \int_\Gamma grad id grad v ds
{
    const Uint idx_f= f.RowIdx->GetIdx();
    IdxT Numb[10];
    
// std::ofstream fil("surf.off");
// fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    LocalP1CL<Point3DCL> GradRef[10], Grad[10];
    P2DiscCL::GetGradientsOnRef( GradRef);

    SMatrixCL<3,3> T;
    double det;
    InterfaceTriangleCL triangle;

    Quad5_2DCL<Point3DCL> e[3];
    for (int i= 0; i<3; ++i)
        e[i]= std_basis<3>( i + 1);

    DROPS_FOR_TRIANG_CONST_TETRA( MG, /*default level*/-1, it)
    {
    	triangle.Init( *it, SmPhi, lsetbnd);

        for (int v= 0; v < 10; ++v)
        { // collect data on all DoF
            const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
            Numb[v]= unk.Exist( idx_f) ? unk( idx_f) : NoIdx;
        }
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T);

        for (int ch= 0; ch < 8; ++ch)
        {
        	triangle.ComputeForChild( ch);
            for (int t= 0; t < triangle.GetNumTriangles(); ++t)
                SF_ImprovedLaplBeltramiOnTriangle( *it, &triangle.GetBary( t),
                		triangle, Grad,  Numb, e, triangle.GetAbsDet( t), f.Data, sf);
        } // Ende der for-Schleife ueber die Kinder
    }
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


bool LevelsetP2CL::Intersects( const TetraCL& t) const
{
    const Uint idx= Phi.RowIdx->GetIdx();
    double PhiV0= Phi.Data[t.GetVertex(0)->Unknowns(idx)];

    for (int i=1; i<4; ++i)
        if( PhiV0*Phi.Data[t.GetVertex(i)->Unknowns(idx)] <= 0) return true;
    for (int i=0; i<6; ++i)
        if( PhiV0*Phi.Data[t.GetEdge(i)->Unknowns(idx)] <= 0) return true;

    return false;
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
    VecDescCL SmPhi= Phi;
    if (curvDiff_>0)
        SmoothPhi( SmPhi.Data, curvDiff_);

    switch (SF_)
    {
      case SF_LB:
        SF_LaplBeltrami( MG_, SmPhi, BndData_, sf_.GetSigma()(std_basis<3>(0), 0.), f); break;
      case SF_Const:
        SF_ConstForce( MG_, SmPhi, BndData_, sf_.GetSigma()(std_basis<3>(0), 0.), f); break;
      case SF_ImprovedLB:
        SF_ImprovedLaplBeltrami( MG_, SmPhi, BndData_, sf_.GetSigma()(std_basis<3>(0), 0.), f); break;
      case SF_ImprovedLBVar:
         SF_ImprovedLaplBeltrami( MG_, SmPhi, BndData_, f, sf_); break;
      default:
        throw DROPSErrCL("LevelsetP2CL::AccumulateBndIntegral not implemented for this SurfaceForceT");
    }
}

double LevelsetP2CL::GetVolume( double translation, bool fine) const
{
    SMatrixCL<3,3> T;
    InterfaceTetraCL tetra;
    double det, absdet, Volume= 0.;
    LocalP2CL<double> ones( 1.);

    for (MultiGridCL::const_TriangTetraIteratorCL it=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(idx.TriangLevel()), end=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(idx.TriangLevel());
        it!=end; ++it) {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        tetra.Init( *it, Phi, BndData_, translation);
        if (fine)
            for (int ch= 0; ch < 8; ++ch) {
                // compute volume
                tetra.ComputeCutForChild( ch);
                Volume+= tetra.quad( ones, absdet, false);
            }
        else {
            // compute volume
            tetra.ComputeCutForChild( 8);
            Volume+= tetra.quad( ones, absdet, false);
        }
    }
#ifdef _PAR
    Volume = ProcCL::GlobalSum(Volume);
#endif
    return Volume;
}

double LevelsetP2CL::AdjustVolume (double vol, double tol, double surface) const
{
    tol*=vol;

    double v0=GetVolume()-vol;
    if (std::abs(v0)<=tol) return 0;

    double d0=0, d1=v0*(surface ? 1.1/surface : 0.23/std::pow(vol,2./3.));
    // Hinweis: surf(Kugel) = [3/4/pi*vol(Kugel)]^(2/3) * 4pi
    double v1=GetVolume(d1)-vol;
    if (std::abs(v1)<=tol) return d1;

    // Sekantenverfahren fuer Startwert
    while (v1*v0 > 0) // gleiches Vorzeichen
    {
        const double d2=d1-1.2*v1*(d1-d0)/(v1-v0);
        d0=d1; d1=d2; v0=v1; v1=GetVolume(d1)-vol;
        if (std::abs(v1)<=tol) return d1;
    }

    // Anderson-Bjoerk fuer genauen Wert
    while (true)
    {
        const double d2=(v1*d0-v0*d1)/(v1-v0),
                     v2=GetVolume(d2)-vol;
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
    double det, absdet;
    InterfacePatchCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
    maxGradPhi= -1.;
    minGradPhi= 1e99;

    double maxNorm;
    double minNorm;

    DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), it)
    {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
        patch.Init( *it, Phi, BndData_);

        // compute maximal norm of grad Phi
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
        {
            gradPhi+= patch.GetPhi(v)*Grad[v];
        }
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

#ifndef _PAR
void LevelsetRepairCL::pre_refine()
/// do nothing
{
}
#else
void LevelsetRepairCL::pre_refine()
/// Tell parallel multigrid about the location of the DOF
{
    GetPMG().AttachTo( &ls_.Phi, &ls_.GetBndData());
}
#endif

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
#ifdef _PAR
    GetPMG().HandleNewIdx(&ls_.idx, &loc_phi);
#endif
    RepairAfterRefineP2( ls_.GetSolution( phi), loc_phi);
#ifdef _PAR
    GetPMG().CompleteRepair( &loc_phi);
#endif

    phi.Clear( phi.t);
    ls_.DeleteNumbering( phi.RowIdx);
    ls_.idx.swap( loc_lidx);
    phi.SetIdx( &ls_.idx);
    phi.Data= loc_phi.Data;
}

} // end of namespace DROPS


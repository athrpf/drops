//**************************************************************************
// File:    levelset.cpp                                                   *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "levelset/levelset.h"
#include "levelset/fastmarch.h"
#include <fstream>

namespace DROPS
{

inline double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha);
}

void SF_ConstForce( const MultiGridCL& MG, const VecDescCL& SmPhi, double sigma, VecDescCL& f)
// computes the integral
//         sigma \int_\Gamma v n ds
// used by levelset/prJump.cpp for testing FE pressure spaces.
{
    const Uint idx_f= f.RowIdx->GetIdx();
    IdxT Numb[10];

//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det;
    InterfacePatchCL patch;
    P2DiscCL::GetGradientsOnRef( GradRef);

    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    for (MultiGridCL::const_TriangTetraIteratorCL it=MG.GetTriangTetraBegin(), end=MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
        patch.Init( *it, SmPhi);

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
            Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
        }

        for (int ch=0; ch<8; ++ch)
        {
            if (!patch.ComputeForChild(ch)) // no patch for this child
                continue;

//patch.WriteGeom( fil);

            BaryCoordCL BaryPQR, BarySQR;
            for (int i=0; i<3; ++i)
            {
                // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
                BaryPQR+= patch.GetBary(i);
                BarySQR+= patch.GetBary(i+1);
            }
            BaryPQR/= 3; BarySQR/= 3;

            Point3DCL n; // senkrecht auf PQ, PR, nach aussen gerichtet...
            cross_product( n, patch.GetPoint(1)-patch.GetPoint(0), patch.GetPoint(2)-patch.GetPoint(0));
            n/= n.norm();

            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            const int find_sign= patch.GetNumSign( 1) ? 1 : -1;
            Point3DCL pos_dir;
            for (Uint i=0; i<4; ++i) // compute vector with positive direction
            {
                const Uint vert= data.Vertices[i];
                if (patch.GetSign(vert)==find_sign)
                {
                    const Point3DCL signedPoint= vert<4 ? it->GetVertex(vert)->GetCoord() : GetBaryCenter( *it->GetEdge(vert-4));
                    pos_dir= signedPoint - patch.GetPoint(0);
                    if (find_sign == -1) pos_dir= -pos_dir;
                    break;
                }
            }
            if (inner_prod( n, pos_dir) < 0) n= -n;

            double val_hat[4];
            for (int v=0; v<10; ++v)
            {
                if (Numb[v]==NoIdx) continue;

                for (Uint k=0; k<patch.GetNumPoints(); ++k)
                    // Werte der Hutfunktion in P,Q,R,S
                    val_hat[k]= FE_P2CL::H(v,patch.GetBary(k));

                double v_Bary= FE_P2CL::H(v,BaryPQR),
                    sum_v= 0;
                for (int k=0; k<3; ++k)
                     sum_v+= val_hat[k];

                if (patch.IsQuadrilateral())
                {
                    double sum_vSQR= 0;
                    for (int k=1; k<4; ++k)
                        sum_vSQR+= val_hat[k];
                    sum_v+= patch.GetAreaFrac() * sum_vSQR;
                    v_Bary+= patch.GetAreaFrac() * FE_P2CL::H(v,BarySQR);
                }

                // Quadraturformel auf Dreieck, exakt bis zum Grad 2
                const double int_v= (1./12)*sum_v + 0.75 * v_Bary;
                const double C= sigma*patch.GetFuncDet()/2;
                for (int i=0; i<3; ++i)
                    f.Data[Numb[v]+i]-= C * int_v*n[i];
            }
        } // Ende der for-Schleife ueber die Kinder
    }
//fil << "}\n";
}

void SF_LaplBeltrami( const MultiGridCL& MG, const VecDescCL& SmPhi, double sigma, VecDescCL& f)
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
    InterfacePatchCL patch;

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

        patch.Init( *it, SmPhi);

        for (int ch=0; ch<8; ++ch)
        {
            if (!patch.ComputeForChild(ch)) // no patch for this child
                continue;

//patch.WriteGeom( fil);

            BaryCoordCL BaryPQR, BarySQR;
            for (int i=0; i<3; ++i)
            {
                // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
                BaryPQR+= patch.GetBary(i);
                BarySQR+= patch.GetBary(i+1);
            }

            const double C= patch.GetFuncDet()/6*sigma;  // 1/6 for quad. rule

            for (int v=0; v<10; ++v)
            {
                if (Numb[v]==NoIdx) continue;

                LocalP1CL<Point3DCL> gradv; // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
                for (int node=0; node<4; ++node)
                    gradv[node]= Grad[v][node];

                Point3DCL gr= gradv( BaryPQR); // gr= grad v(P) + grad v(Q) + grad v(R)

                if (patch.IsQuadrilateral())
                    gr+= patch.GetAreaFrac() * gradv( BarySQR);
                // nun gilt:
                // gr = [grad v(P)+...] + (a+b-1)[grad v(S)+...]

                for (int i=0; i<3; ++i)
                {
                    const double val= inner_prod( gr, patch.GetGradId(i));
                    f.Data[Numb[v]+i]-= C *val;
                }
            }
        } // Ende der for-Schleife ueber die Kinder
    }
//fil << "}\n";
}

void SF_ImprovedLaplBeltrami( const MultiGridCL& MG, const VecDescCL& SmPhi, double sigma, VecDescCL& f)
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
  InterfacePatchCL patch;

  P2DiscCL::GetGradientsOnRef( GradRef);

  for (MultiGridCL::const_TriangTetraIteratorCL it=MG.GetTriangTetraBegin(), end=MG.GetTriangTetraEnd();
       it!=end; ++it)
  {
    GetTrafoTr( T, det, *it);
    P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
    LocalP1CL<Point3DCL> n;

    patch.Init( *it, SmPhi);
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
      const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
      Numb[v]= unk.Exist(idx_f) ? unk(idx_f) : NoIdx;
      for (int k=0; k<4; ++k)
        n[k]+= patch.GetPhi(v)*Grad[v][k];
    }

    for (int ch=0; ch<8; ++ch)
    {
      if (!patch.ComputeForChild(ch)) // no patch for this child
        continue;

//patch.WriteGeom( fil);
      BaryCoordCL BaryPQR, BarySQR;
      for (int i=0; i<3; ++i)
      {
                // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
        BaryPQR+= patch.GetBary(i);
        BarySQR+= patch.GetBary(i+1);
      }
      BaryPQR/= 3.;    BarySQR/= 3.;

      typedef SArrayCL<Point3DCL,3> ProjT;
      GridFunctionCL<ProjT> GradId( ProjT(), 6);  // values in P, Q, R, S, BaryPQR, BarySQR
      for (int p=0; p<6; ++p)
      {
        Point3DCL np= n( p<4 ? patch.GetBary(p) : p==4 ? BaryPQR : BarySQR);
        if (np.norm()>1e-8) np/= np.norm();
        for (int i=0; i<3; ++i)
          GradId[p][i]= patch.ApplyProj( std_basis<3>(i+1) - np[i]*np);
//                     GradId[p][i]= std_basis<3>(i+1) - np[i]*np;
      }

      const double C= patch.GetFuncDet()*sigma/2.;

      for (int v=0; v<10; ++v)
      {
        if (Numb[v]==NoIdx) continue;

        LocalP1CL<Point3DCL> gradv; // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
        for (int node=0; node<4; ++node)
          gradv[node]= Grad[v][node];

        for (int i=0; i<3; ++i)
        {
          double intSum= 0; // sum of the integrand in PQR, SQR
          for (int k=0; k<3; ++k)
          {
            intSum+= inner_prod( GradId[k][i], gradv(patch.GetBary(k)));
            if (patch.IsQuadrilateral())
              intSum+= patch.GetAreaFrac() * inner_prod( GradId[k+1][i], gradv(patch.GetBary(k+1)));
          }
          double intBary= inner_prod( GradId[4][i], gradv(BaryPQR));
          if (patch.IsQuadrilateral())
            intBary+= patch.GetAreaFrac() * inner_prod( GradId[5][i], gradv(BarySQR));
          f.Data[Numb[v]+i]-= C *(intSum/12. + 0.75*intBary);
        }
      }
    } // Ende der for-Schleife ueber die Kinder
  }
//fil << "}\n";
}

void SF_ImprovedLaplBeltramiOnTriangle( const TetraCL& t, const BaryCoordCL * const p,
                                        const InterfacePatchCL&  patch, const LocalP1CL<Point3DCL> Grad_f[10], const IdxT Numb[10],
                                        instat_scalar_fun_ptr sigma, const Quad5_2DCL<Point3DCL> e[3],
                                        double det, VectorCL& f)
{
    static Quad5_2DCL<Point3DCL> Grad[10]; // Gradients of the P2-basis-functions
    Quad5_2DCL<Point3DCL> n;

    for (int v=0; v<10; ++v)
    {
        Grad[v].assign( Grad_f[v], p);
        n+= patch.GetPhi(v)*Grad[v];
    }
    for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();

    Quad5_2DCL<> qsigma( t, p, sigma),  // surface tension
                 q1;                    // Term 1

    Quad5_2DCL<Point3DCL> qPhPhte,      // Common term in Term 1 and Term 2
                          qsigmaPhPhte; // for Term 1

    for (int i= 0; i < 3; ++i)
    {
        qPhPhte= (e[i] - dot(e[i],n)*n);
        qPhPhte.apply( patch, &InterfacePatchCL::ApplyProj);
        qsigmaPhPhte= qsigma*qPhPhte;
        for (int v= 0; v < 10; ++v)
        {
            q1= dot (qsigmaPhPhte, Grad[v]);
            f[Numb[v]+i]-= q1.quad( det);
        }
    }
}

void SF_ImprovedLaplBeltramiOnTriangle( const TetraCL& t, const BaryCoordCL * const p,
    const InterfacePatchCL&  patch, const LocalP1CL<Point3DCL> Grad_f[10], const IdxT Numb[10],
    instat_scalar_fun_ptr sigma, instat_vector_fun_ptr grad_sigma, const Quad5_2DCL<Point3DCL> e[3],
    double det, VectorCL& f)
{
    if (grad_sigma == 0)
    {
        SF_ImprovedLaplBeltramiOnTriangle( t, p, patch, Grad_f, Numb, sigma, e, det, f);
        return;
    }
    static Quad5_2DCL<>          p2[10];   // P2-Hat-Functions...
    static Quad5_2DCL<Point3DCL> Grad[10]; // and their gradients
    Quad5_2DCL<Point3DCL> n;
    P2DiscCL::GetP2Basis( p2, p);
    for (int v=0; v<10; ++v)
    {
        Grad[v].assign( Grad_f[v], p);
        n+= patch.GetPhi(v)*Grad[v];
    }
    for (int i =0; i<Quad5_2DDataCL::NumNodesC; i++) if (n[i].norm()>1e-8) n[i]/= n[i].norm();

    Quad5_2DCL<> qsigma( t, p, sigma), // surface tension
                 q1,                   // Term 1
                 q2_3;                 // Term 2 minus Term 3
    Quad5_2DCL<Point3DCL> qPhPhte,                         // Common term in Term 1 and Term 2
                          qsigmaPhPhte,                    // for Term 1
                          qgradsigma( t, p, grad_sigma),   // for Term 2
                          qPhgradsigma( qgradsigma);       // for Term 3
    Quad5_2DCL<>          qgradsigmaPhPhte_m_Phgradsigmae; // for Term 2 and 3
    qPhgradsigma.apply( patch, &InterfacePatchCL::ApplyProj);

    for (int i= 0; i < 3; ++i)
    {
        qPhPhte= (e[i] - dot(e[i],n)*n);
        qPhPhte.apply( patch, &InterfacePatchCL::ApplyProj);

        qsigmaPhPhte= qsigma*qPhPhte;
        qgradsigmaPhPhte_m_Phgradsigmae= dot( qPhPhte, qgradsigma) - dot( qPhgradsigma, e[i]);
        for (int v= 0; v < 10; ++v)
        {
            q1= dot (qsigmaPhPhte, Grad[v]);
            q2_3= p2[v]*qgradsigmaPhPhte_m_Phgradsigmae;
            f[Numb[v]+i]-= q1.quad( det) + q2_3.quad( det);
        }
    }
}

void SF_ImprovedLaplBeltrami( const MultiGridCL& MG, const VecDescCL& SmPhi,
    instat_scalar_fun_ptr sigma, instat_vector_fun_ptr grad_sigma, VecDescCL& f)
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
    InterfacePatchCL patch;

    Quad5_2DCL<Point3DCL> e[3];
    for (int i= 0; i<3; ++i)
        e[i]= std_basis<3>( i + 1);

    DROPS_FOR_TRIANG_CONST_TETRA( MG, /*default level*/, it)
    {
        patch.Init( *it, SmPhi);

        for (int v= 0; v < 10; ++v)
        { // collect data on all DoF
            const UnknownHandleCL& unk= v<4 ? it->GetVertex(v)->Unknowns : it->GetEdge(v-4)->Unknowns;
            Numb[v]= unk.Exist( idx_f) ? unk( idx_f) : NoIdx;
        }
        GetTrafoTr( T, det, *it);
        P2DiscCL::GetGradients( Grad, GradRef, T);

        for (int ch= 0; ch < 8; ++ch)
        {
            patch.ComputeForChild( ch);
            for (int t= 0; t < patch.GetNumTriangles(); ++t)
                SF_ImprovedLaplBeltramiOnTriangle( *it, &patch.GetBary( t),
                    patch, Grad,  Numb, sigma, grad_sigma, e, patch.GetFuncDet( t), f.Data);
        } // Ende der for-Schleife ueber die Kinder
    }
}

void MarkInterface ( scalar_fun_ptr DistFct, double width, MultiGridCL& mg)
{
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/, it)
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
    DROPS_FOR_TRIANG_TETRA( mg, /*default-level*/, it)
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
        Phi.Data[it->Unknowns(idx)]= phi0( it->GetCoord());
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= MG_.GetTriangEdgeBegin(lvl),
        end= MG_.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( GetBaryCenter( *it));
    }
}

void LevelsetP2CL::Reparam( Uint steps, double dt)
// Reparametrization of the levelset function Phi
{
    VectorCL Psi= Phi.Data, b;
    MatrixCL L, R, M;

    for (Uint i=0; i<steps; ++i)
    {
        SetupReparamSystem( M, R, Psi, b);
        L.LinComb( 1./dt, M, theta_, R);

        b+= (1./dt)*(M*Psi) - (1.-theta_) * (R*Psi);
        gm_.Solve( L, Psi, b);
        std::cout << "Reparam: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
    }

    Phi.Data= Psi;
}

double func_abs( double x) { return std::abs(x); }

void LevelsetP2CL::SetupReparamSystem( MatrixCL& M_, MatrixCL& R_, const VectorCL& Psi, VectorCL& b) const
// M, R, b describe the following terms used for reparametrization:
// b_i  = ( S(Phi0),           v_i     )
// M_ij = ( v_j,               v_i     )
// R_ij = ( w(Psi) grad v_j,   v_i     )
//      + (|S(Phi0)| grad v_j, grad v_i) * diff
// where v_i, v_j denote the ansatz functions
// and w(Psi) = sign(Phi0) * grad Psi / |grad Psi| the scaled gradient of Psi
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.GetLevel();

    SparseMatBuilderCL<double> R(&R_, num_unks, num_unks);
    SparseMatBuilderCL<double> M(&M_, num_unks, num_unks);
    b.resize( 0);
    b.resize( num_unks);
    std::cerr << "entering SetupReparamSystem: " << num_unks << " levelset unknowns. ";

    Quad2CL<double>    Sign_Phi;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], w_loc;
    SMatrixCL<3,3>     T;
    P2DiscCL::GetGradientsOnRef( GradRef);

    IdxT         Numb[10];
    SVectorCL<3> grad_Psi[4];
    const_DiscSolCL    phi= GetSolution();
    double det, absdet;
    const double alpha= 0.1;  // for smoothing of signum fct

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);

        GetLocalNumbP2NoBnd( Numb, *sit, *Phi.RowIdx);

        // init Sign_Phi, w_loc
        for (int i=0; i<4; ++i)
        {
            Sign_Phi[i]= SmoothedSign( phi.val( *sit->GetVertex(i)), alpha);
            grad_Psi[i]= Point3DCL();  // init with zero
            for (int l=0; l<10; ++l)
                grad_Psi[i]+= Psi[ Numb[l]] * Grad[l][i];
            w_loc[i]= (Sign_Phi[i]/grad_Psi[i].norm() )*grad_Psi[i];
        }
        // values in barycenter
        Point3DCL gr= 0.25*(grad_Psi[0]+grad_Psi[1]+grad_Psi[2]+grad_Psi[3]);
        Sign_Phi[4]= SmoothedSign( phi.val( *sit, 0.25, 0.25, 0.25), alpha);
        w_loc[4]=    Sign_Phi[4]*gr/gr.norm();

        Quad2CL<> w_Grad[10];
        for (int i=0; i<10; ++i)
            w_Grad[i]= dot(w_loc, Grad[i]);

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            // b_i  = ( S(Phi0),         v_i + SD * w(Psi) grad v_i )
            b[ Numb[i]]+= Sign_Phi.quadP2( i, absdet);
            for(int j=0; j<10; ++j)
            {
                // M_ij = ( v_j, v_i )
                M( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j)*absdet;
                // R_ij = ( w(Psi) grad v_j, v_i + SD * w(Psi) grad v_i )
                R( Numb[i], Numb[j])+= w_Grad[j].quadP2(i, absdet)
                    + diff_*Quad2CL<>(dot( Grad[j]*Sign_Phi.apply( func_abs), Grad[i])).quad( absdet);
            }
        }
    }
    M.Build();
    R.Build();
    std::cerr << M_.num_nonzeros() << " nonzeros in M, "
              << R_.num_nonzeros() << " nonzeros in R!" << std::endl;
}

void LevelsetP2CL::SetTimeStep( double dt, double theta)
{
    dt_= dt;
    if (theta >= 0) theta_= theta;

    L_.LinComb( 1./dt_, E, theta_, H);
}

void LevelsetP2CL::ComputeRhs( VectorCL& rhs) const
{
std::cerr << "ComputeRhs-dt_: " << dt_ << std::endl;
    rhs= (1./dt_)*Phi.Data;
    if (theta_ != 1.) {
        GMResSolverCL<SSORPcCL> gm( gm_);
        VectorCL tmp( rhs.size());
        gm.Solve( E, tmp, VectorCL( H*Phi.Data));
        std::cerr << "ComputeRhs: res = " << gm.GetResid() << ", iter = " << gm.GetIter() << std::endl;
        rhs-= (1. - theta_)*tmp;
    }
}

void LevelsetP2CL::DoLinStep( const VectorCL& rhs)
{
    gm_.Solve( L_, Phi.Data, rhs);
    std::cerr << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() <<std::endl;
}

void LevelsetP2CL::DoStep( const VectorCL& rhs)
{
    gm_.Solve( L_, Phi.Data, VectorCL( E*rhs));
    std::cerr << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() <<std::endl;
}

void LevelsetP2CL::DoStep()
{
    VectorCL rhs( Phi.Data.size());
    ComputeRhs( rhs);
    DoStep( rhs);
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


void LevelsetP2CL::ReparamFastMarching( bool ModifyZero, bool Periodic, bool OnlyZeroLvl)
/** \param ModifyZero: If true, the zero level is moved inside the elements intersecting the interface. If false, the zero level is kept fixed.
 *  \param OnlyZeroLvl: If true, only the first step of the algorithm is performed, i.e. the reparametrization only takes place locally at the interface.
 *  \param Periodic: If true, a special variant of the algorithm for periodic boundaries is used.
 */
{
    FastMarchCL fm( MG_, Phi);

    if (OnlyZeroLvl)
    {
        if (Periodic)
            fm.InitZeroPer( Bnd_, ModifyZero);
        else
            fm.InitZero( ModifyZero);
        fm.RestoreSigns();
    }
    else if (Periodic)
        fm.ReparamPer( Bnd_, ModifyZero);
    else
        fm.Reparam( ModifyZero);
}


void LevelsetP2CL::AccumulateBndIntegral( VecDescCL& f) const
{
    VecDescCL SmPhi= Phi;
    if (curvDiff_>0)
        SmoothPhi( SmPhi.Data, curvDiff_);
    switch (SF_)
    {
      case SF_LB:
        SF_LaplBeltrami( MG_, SmPhi, sigma(std_basis<3>(0), 0.), f); break;
      case SF_Const:
        SF_ConstForce( MG_, SmPhi, sigma(std_basis<3>(0), 0.), f); break;
      case SF_ImprovedLB:
        SF_ImprovedLaplBeltrami( MG_, SmPhi, sigma(std_basis<3>(0), 0.), f); break;
      case SF_ImprovedLBVar:
        SF_ImprovedLaplBeltrami( MG_, SmPhi, sigma, grad_sigma, f); break;
      default:
        throw DROPSErrCL("LevelsetP2CL::AccumulateBndIntegral not implemented for this SurfaceForceT");
    }
}

double LevelsetP2CL::GetVolume( double translation) const
{
    SMatrixCL<3,3> T;
    InterfacePatchCL patch;
    double det, absdet, Volume= 0.;
    LocalP2CL<double> ones( 1.); 
    
    for (MultiGridCL::const_TriangTetraIteratorCL it=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(), end=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd();
        it!=end; ++it) {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        patch.Init( *it, Phi, translation);
        for (int ch= 0; ch < 8; ++ch) {
            // compute volume
            patch.ComputeCutForChild( ch);
            Volume+= patch.quad( ones, absdet, false);
        }
    }
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
    std::cerr << "Smoothing for curvature calculation... ";
    MatrixCL M, A, C;
    SetupSmoothSystem( M, A);
    C.LinComb( 1, M, diff, A);
    PCG_SsorCL pcg( pc_, gm_.GetMaxIter(), gm_.GetTol());
    pcg.Solve( C, SmPhi, M*Phi.Data);
    std::cerr << "||SmPhi - Phi||_oo = " << supnorm( SmPhi-Phi.Data) << std::endl;
}

void LevelsetP2CL::SetupSmoothSystem( MatrixCL& M, MatrixCL& A) const
// used for smoothing of Phi before computing curvature term
//
// M = mass matrix for P2 elements
// A = stiffness matrix for P2 elements
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.GetLevel();

    SparseMatBuilderCL<double> Mb(&M, num_unks, num_unks);
    SparseMatBuilderCL<double> Ab(&A, num_unks, num_unks);

    Quad2CL<Point3DCL> Grad[10], GradRef[10], w_loc;
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

//*****************************************************************************
//                               LevelsetRepairCL
//*****************************************************************************
inline void
LevelsetRepairCL::post_refine ()
{
    VecDescCL loc_phi;
    IdxDescCL loc_lidx( P2_FE);
    VecDescCL& phi= ls_.Phi;
    match_fun match= ls_.GetMG().GetBnd().GetMatchFun();

    ls_.CreateNumbering( ls_.GetMG().GetLastLevel(), &loc_lidx, match);
    loc_phi.SetIdx( &loc_lidx);
    RepairAfterRefineP2( ls_.GetSolution( phi), loc_phi);

    phi.Clear();
    ls_.DeleteNumbering( phi.RowIdx);
    ls_.idx.swap( loc_lidx);
    phi.SetIdx( &ls_.idx);
    phi.Data= loc_phi.Data;
}

//*****************************************************************************
//                               InterfacePatchCL
//*****************************************************************************

const double InterfacePatchCL::approxZero_= 2.*std::numeric_limits<double>::epsilon();
const bool   InterfacePatchCL::LinearEdgeIntersection;
BaryCoordCL  InterfacePatchCL::AllEdgeBaryCenter_[10][10];
BaryCoordCL  InterfacePatchCL::BaryDoF_[10];

InterfacePatchCL::InterfacePatchCL()
  : RegRef_( GetRefRule( RegRefRuleC)), intersec_(0), ch_(-1)
{
    BaryDoF_[0][0]= BaryDoF_[1][1]= BaryDoF_[2][2]= BaryDoF_[3][3]= 1.;
    for (int edge=0; edge<6; ++edge)
        BaryDoF_[edge+4]= 0.5*(BaryDoF_[VertOfEdge(edge,0)] + BaryDoF_[VertOfEdge(edge,1)]);
    for (int i= 0; i < 10; ++i)
        for (int j= 0; j < 10; ++j)
            AllEdgeBaryCenter_[i][j]= BaryCenter( BaryDoF_[i], BaryDoF_[j]);
}

void InterfacePatchCL::Init( const TetraCL& t, const VecDescCL& ls, double translation)
{
    const Uint idx_ls= ls.RowIdx->GetIdx();
    ch_= -1;
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        Coord_[v]= v<4 ? t.GetVertex(v)->GetCoord() : GetBaryCenter( *t.GetEdge(v-4));
        PhiLoc_[v]= ls.Data[v<4 ? t.GetVertex(v)->Unknowns(idx_ls) : t.GetEdge(v-4)->Unknowns(idx_ls)] + translation;
        sign_[v]= Sign(PhiLoc_[v]);
    }
}

void InterfacePatchCL::Init( const TetraCL& t, const LocalP2CL<double>& ls, double translation)
{
    ch_= -1;
    PhiLoc_= ls + translation;
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        Coord_[v]= v<4 ? t.GetVertex(v)->GetCoord() : GetBaryCenter( *t.GetEdge(v-4));
        sign_[v]= Sign(PhiLoc_[v]);
    }
}

bool InterfacePatchCL::ComputeForChild( Uint ch)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
    ch_= ch;
    num_sign_[0]= num_sign_[1]= num_sign_[2]= 0;
    for (int vert= 0; vert<4; ++vert)
        ++num_sign_[ sign_[data.Vertices[vert]] + 1];

    intersec_= 0;
    if (num_sign_[0]*num_sign_[2]==0 && num_sign_[1]<3) {// no change of sign on child
        numchildtriangles_= 0;
        return false;
    }
    if (num_sign_[1]==4)
    {
        std::cerr << "WARNING: InterfacePatchCL: found 3-dim. zero level set, grid is too coarse!" << std::endl;
        numchildtriangles_= 0;
        return false;
    }

    // erst werden die Nullknoten in PQRS gespeichert...
    for (int vert= 0; vert<4; ++vert)
    {
        const int v= data.Vertices[vert];
        if (sign_[v]==0)
        {
            Bary_[intersec_]= BaryDoF_[v];
            PQRS_[intersec_++]= Coord_[v];
        }
    }
    // ...dann die echten Schnittpunkte auf den Kanten mit Vorzeichenwechsel
    for (int edge= 0; edge<6; ++edge)
    {
        const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                  v1= data.Vertices[ VertOfEdge( edge, 1)];
        if (sign_[v0]*sign_[v1]<0) // different sign -> 0-level intersects this edge
        {
            const double lambda= EdgeIntersection( v0,  v1, PhiLoc_);
            Bary_[intersec_]= (1-lambda)*BaryDoF_[v0] + lambda * BaryDoF_[v1];
            // bary-coords of tetra, not of subtetra!
            PQRS_[intersec_++]= (1-lambda) * Coord_[v0] + lambda * Coord_[v1];
        }
    }
    if (intersec_<3) { // Nullstellenmenge vom Mass 0!
        numchildtriangles_= 0;
        return false;
    }

    SMatrixCL<3,2> A;    // A = [ Q-P | R-P ]
    A(0,0)= PQRS_[1][0]-PQRS_[0][0];    A(0,1)= PQRS_[2][0]-PQRS_[0][0];
    A(1,0)= PQRS_[1][1]-PQRS_[0][1];    A(1,1)= PQRS_[2][1]-PQRS_[0][1];
    A(2,0)= PQRS_[1][2]-PQRS_[0][2];    A(2,1)= PQRS_[2][2]-PQRS_[0][2];
    SMatrixCL<2,2> ATA;
    ATA(0,0)=           A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
    ATA(0,1)= ATA(1,0)= A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
    ATA(1,1)=           A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
    const double detATA= ATA(0,0)*ATA(1,1) - ATA(1,0)*ATA(1,0);
    sqrtDetATA_= std::sqrt( detATA);

    Point2DCL AT_i, tmp;
    for (int i=0; i<3; ++i)
    {
        // berechne B = A * (ATA)^-1 * AT
        AT_i[0]= A(i,0); AT_i[1]= A(i,1);
        Solve2x2( detATA, ATA, tmp, AT_i);
        B_[i]= A*tmp;
    }

    if (intersec_==4) // 4 intersections --> a+b != 1
    { // berechne a, b
        // Loese (Q-P)a + (R-P)b = S-P  --> lin. AGP, loese ATA * [a,b]T = AT(S-P)
        Point3DCL PS= PQRS_[3] - PQRS_[0];
        tmp[0]= A(0,0)*PS[0] + A(1,0)*PS[1] + A(2,0)*PS[2];
        tmp[1]= A(0,1)*PS[0] + A(1,1)*PS[1] + A(2,1)*PS[2];
        Solve2x2( detATA, ATA, ab_, tmp);
        //if (ab_[0]<0 || ab_[1]<0)
        //    std::cerr<<"LevelsetP2CL::AccumulateBndIntegral: a or b negative"<<std::endl;
        // a,b>=0 muss erfuellt sein, da wegen edge+oppEdge==5 die Punkte P und S sich automatisch gegenueber liegen muessten...
        numchildtriangles_= 2;
    }
    else
        numchildtriangles_= 1;

    if (EqualToFace()) // interface is shared by two tetras
        sqrtDetATA_/= 2;

    return true; // computed patch of child;
}


Point3DCL InterfacePatchCL::GetNormal() const
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch_]);
    SMatrixCL<3,4> p1grad;
    double det; // dummy
    Point3DCL pt[4];
    SVectorCL<4> ls;
    for (int v= 0; v < 4; ++v) {
        pt[v]= Coord_ [data.Vertices[v]];
        ls[v]= PhiLoc_[data.Vertices[v]];
    }
    P1DiscCL::GetGradients( p1grad, det, pt);
    const Point3DCL n( p1grad*ls);
    return n/n.norm();
}


bool InterfacePatchCL::ComputeCutForChild( Uint ch)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
    ch_= ch;
    num_sign_[0]= num_sign_[1]= num_sign_[2]= 0;
    for (int vert= 0; vert<4; ++vert)
        ++num_sign_[ sign_[data.Vertices[vert]] + 1];

    intersec_= innersec_= 0;
    if (num_sign_[0]*num_sign_[2]==0 && num_sign_[1]<3) // no change of sign on child and no patch on a face
        return false;
    if (num_sign_[1]==4)
    {
        std::cerr << "WARNING: InterfacePatchCL: found 3-dim. zero level set, grid is too coarse!" << std::endl;
        return false;
    }

    // erst werden die Nullknoten in PQRS gespeichert...
    for (int vert= 0; vert<4; ++vert)
    {
        const int v= data.Vertices[vert];
        if (sign_[v]==0)
        {
            Bary_[intersec_]= BaryDoF_[v];
            Edge_[intersec_++]= -1;
        }
    }
    // ...dann die echten Schnittpunkte auf den Kanten mit Vorzeichenwechsel
    for (int edge= 0; edge<6; ++edge)
    {
        const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                  v1= data.Vertices[ VertOfEdge( edge, 1)];
        if (sign_[v0]*sign_[v1]<0) // different sign -> 0-level intersects this edge
        {
            const double lambda= EdgeIntersection( v0,  v1, PhiLoc_);
            Bary_[intersec_]= (1-lambda)*BaryDoF_[v0] + lambda * BaryDoF_[v1];
            Edge_[intersec_++]= edge;
            innersec_++;
        }
    }
    if (intersec_<3) return false; // Nullstellenmenge vom Mass 0!

    return true; // computed cut of child;
}

void InterfacePatchCL::DebugInfo( std::ostream& os, bool InfoForChild) const
{
    if (InfoForChild)
    {
        const ChildDataCL data= GetChildData( RegRef_.Children[ch_]);
        os << "Patch on child " << ch_ << " with " << GetNumPoints() << " intersections:\n";
        for (Uint i=0; i<GetNumPoints(); ++i)
            os << "\tP" << i << " = " << GetPoint(i) << '\n';
        os << "Signs of verts: ";
        for (int i=0; i<4; ++i)
            os << GetSign(data.Vertices[i]) << " ";
        os << std::endl;
    }
    else
    {
        os << "Signs: ";
        for (int i=0; i<10; ++i)
        {
            os << GetSign(i) << " ";
        }
        os << std::endl;
    }
}

void InterfacePatchCL::WriteGeom( std::ostream& os) const
{
    os << "geom {OFF " << intersec_ << " 1 0\n";
    for (int i=0; i<intersec_; ++i)
    {
        for (int j=0; j<3; ++j)
            os << PQRS_[i][j] << ' ';
        os << '\n';
    }
    if (IsQuadrilateral())
        os << "4 0 1 3 2";
    else
        os << "3 0 1 2";
    os << "\n}\n";
}

void InterfacePatchCL::InsertSubTetra(SubTetraT& BaryCoords, bool pos)
{
    if (pos)
        posTetras.push_back(BaryCoords);
    else
        negTetras.push_back(BaryCoords);
}

void InterfacePatchCL::ComputeSubTets()
{
    posTetras.clear();
    negTetras.clear();

    ChildDataCL data;
    SubTetraT BaryCoords;

    if (!Intersects())
    {
        for (Uint i=0; i<4; ++i)
            BaryCoords[i] = BaryDoF_[i];
        InsertSubTetra( BaryCoords, sign_[0]==1);
        return;
    }

    // Schleife ueber die Kinder
    for (Uint ch=0; ch < 8 && Intersects(); ++ch)
    {
        //std::cerr << "Kind " << ch << "\n \n";
        data  = GetChildData (RegRef_.Children[ch]);

        //cuts = Child + empty set
        if (!ComputeCutForChild(ch) || (intersec_==3 && innersec_==0))
        {
            //std::cerr << "cuts = Child + empty set: " << num_sign_[2] <<"\n";
            for (Uint i=0; i<4; ++i)
                BaryCoords[i]=BaryDoF_[data.Vertices[i]];
            InsertSubTetra( BaryCoords, num_sign_[2]>0);
        }
        if (intersec_==3)
            switch (innersec_)
            {
                case 1 : // cuts = Tetra + Tetra
                {
                    //std::cerr << "cuts = Tetra + Tetra"<<"\n";
                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];
                    for (int i=0; i<4; ++i)
                    {
                        if (sign_[data.Vertices[i]]==0) continue;
                        BaryCoords[3]= BaryDoF_[data.Vertices[i]];
                        InsertSubTetra( BaryCoords, sign_[data.Vertices[i]]==1);
                    }
                } break;
                case 2 : // cuts = Tetra + Pyramide mit 4eckiger Grundflaeche
                {
                    // tetra
                    //std::cerr << "cuts = Tetra + Pyramide"<<"\n";
                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];

                    int vertB= -1;
                    const int signB= num_sign_[0]==1 ? -1 : 1; // sign of vert B occurs only once
                    for (int i=0; vertB==-1 && i<4; ++i)
                        if (sign_[data.Vertices[i]]==signB) vertB= i;
                    BaryCoords[3]= BaryDoF_[data.Vertices[vertB]];
                    InsertSubTetra( BaryCoords, signB==1);

                    // pyramid = 2 tetras: ACDP + APQD
                    //                                     connectivity:     P-------Q
                    //                                                       | \   / |
                    //                                                       |   A   |
                    //                                                       | /   \ |
                    //                                                       C-------D
                    //A,C,D Ecken des Kindtetraeders
                    int z=0;
                    for (int i=0; i<4; ++i)
                        if (i!=vertB) BaryCoords[z++]= BaryDoF_[data.Vertices[i]];
                    BaryCoords[3]= Bary_[1];
                    InsertSubTetra( BaryCoords, signB!=1);

                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];
                    int vertD=-1;
                    for (int i=0; vertD==-1 && i<4; ++i)
                    {
                        if (sign_[data.Vertices[i]]==-signB && 
                            (i==VertOfEdge(Edge_[2],0) ||
                             i==VertOfEdge(Edge_[2],1))  )
                             vertD= i;
                    }
                    BaryCoords[3]= BaryDoF_[data.Vertices[vertD]];
                    InsertSubTetra( BaryCoords, signB!=1);
                } break;
                case 3 : // cuts = Tetra + Tetra-Stumpf
                {
                    //std::cerr << "cuts = Tetra + Tetra-Stumpf\n";
                    int vertA= -1;  // cut-Tetra = APQR
                    const int signA= num_sign_[0]==1 ? -1 : 1; // sign of vert A occurs only once
                    for (int i=0; vertA==-1 && i<4; ++i)
                        if (sign_[data.Vertices[i]]==signA) vertA= i;
                    for (int i=0; i<3; ++i)
                        BaryCoords[i]= Bary_[i];
                    BaryCoords[3]= BaryDoF_[data.Vertices[vertA]];

                    InsertSubTetra( BaryCoords, signA==1);
                    //                                     connectivity:     R---------D
                    //                                                      /|        /|
                    //                                                     Q-|-------C |
                    //                                                      \|        \|
                    //                                                       P---------B
                    //B,C,D Ecken des Kindtetraeders
                    int vertBCD[3];
                    for (Uint i=0; i<3; ++i)
                        vertBCD[i]= (vertA==VertOfEdge(Edge_[i],0)) ? VertOfEdge(Edge_[i],1) : VertOfEdge(Edge_[i],0);

                    // QCPR
                    BaryCoords[0] = Bary_[1];
                    BaryCoords[1] = BaryDoF_[data.Vertices[vertBCD[1]]];
                    BaryCoords[2] = Bary_[0];
                    BaryCoords[3] = Bary_[2];
                    if ( signA!=1 && Edge_[1]!=-1)
                        posTetras.push_back(BaryCoords);
                    else
                        if (Edge_[1]!=-1) negTetras.push_back(BaryCoords);

                    // BCPR
                    BaryCoords[0] = BaryDoF_[data.Vertices[vertBCD[0]]];
                    if ( signA!=1 && Edge_[0]!=-1)
                        posTetras.push_back(BaryCoords);
                    else
                        if (Edge_[0]!=-1) negTetras.push_back(BaryCoords);

                    // BCDR
                    BaryCoords[2] = BaryDoF_[data.Vertices[vertBCD[2]]];
                    if ( signA!=1 && Edge_[2]!=-1)
                        posTetras.push_back(BaryCoords);
                    else
                        if (Edge_[2]!=-1) negTetras.push_back(BaryCoords);
                } break;
            }
        if (intersec_==4)
        {   // cuts = 2x 5-Flaechner mit 6 Knoten. connectivity:     R---------S
            //                                                      /|        /|
            //                                                     A-|-------B |
            //                                                      \|        \|
            //                                                       P---------Q
            //std::cerr << "cuts = 2x 5-Flaechner\n";
            int vertAB[2]; // cut mit VZ==part = ABPQRS

            //erst werden die "negativen" Kinder in die Liste hinzugefuegt, dann die "positiven"
            for (int signAB = -1; signAB<=1; signAB+=2) {
                for (int i=0, k=0; i<4 && k<2; ++i)
                    if (sign_[data.Vertices[i]]==signAB) vertAB[k++]= i;
                // connectivity AP automatisch erfuellt, check for connectivity AR
                const bool AR= vertAB[0]==VertOfEdge(Edge_[2],0) || vertAB[0]==VertOfEdge(Edge_[2],1);
                // Integriere ueber Tetras ABPR, QBPR, QBSR    (bzw. mit vertauschten Rollen von Q/R)
                // ABPR    (bzw. ABPQ)
                BaryCoords[0]= BaryDoF_[data.Vertices[vertAB[0]]];
                BaryCoords[1]= BaryDoF_[data.Vertices[vertAB[1]]];
                BaryCoords[2]= Bary_[0];    BaryCoords[3]= Bary_[AR ? 2 : 1];
                InsertSubTetra( BaryCoords, signAB!=-1);
                // QBPR    (bzw. RBPQ)
                BaryCoords[0]=Bary_[AR ? 1 : 2];
                InsertSubTetra( BaryCoords, signAB!=-1);
                // QBSR    (bzw. RBSQ)
                BaryCoords[2]=Bary_[3];
                InsertSubTetra( BaryCoords, signAB!=-1);
            }
        } //intersec_==4 Ende
    } //Ende der Schleife ueber die Kinder
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

    DROPS_FOR_TRIANG_TETRA( MG_, MG_.GetLastLevel(), it)
    {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
        patch.Init( *it, Phi);

        // compute maximal norm of grad Phi
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
        {
            gradPhi+= patch.GetPhi(v)*Grad[v];
        }
        VectorCL normGrad( 5); 
        for (int v=0; v<5; ++v) // init normGrad
            normGrad[v]= norm( gradPhi[v]);
        const double maxNorm= normGrad.max();
        const double minNorm= normGrad.min();
        if (maxNorm > maxGradPhi) maxGradPhi= maxNorm;
        if (minNorm < minGradPhi) minGradPhi= minNorm;
    }
}

} // end of namespace DROPS


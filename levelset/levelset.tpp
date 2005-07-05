//**************************************************************************
// File:    levelset.tpp                                                   *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "num/discretize.h"
#include "levelset/fastmarch.h"
#include <fstream>

namespace DROPS
{

inline double Sign( double x)
{
    return x<0 ? -1 : x>0 ? 1 : 0;
}

inline double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha);
}

void LevelsetP2CL::DeleteNumbering( IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetAllVertexBegin(level), _MG.GetAllVertexEnd(level) );
    DeleteNumbOnSimplex( idxnum, _MG.GetAllEdgeBegin(level), _MG.GetAllEdgeEnd(level) );
}

void LevelsetP2CL::Init( scalar_fun_ptr phi0)
{
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();
	 
    for (MultiGridCL::TriangVertexIteratorCL it= _MG.GetTriangVertexBegin(lvl),
        end= _MG.GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( it->GetCoord());
    }
    for (MultiGridCL::TriangEdgeIteratorCL it= _MG.GetTriangEdgeBegin(lvl),
        end= _MG.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( GetBaryCenter( *it));
    }
}

/*
inline double QuadP2( double f[10], int i)
{
    double sum, result;
    if (i<4) // hat function on vert
    {
        // Q = sum c[i]*f[i]
        // Gewichte c[i] = 1/420 	fuer vert i
        //                 1/2520	fuer uebrige verts
        //                -1/630        fuer an vert i anliegende edges
        //                -1/420        fuer uebrige drei edges   
        result= f[i]*(1/420.-1./2520.);
        sum= 0.;
        for (int k=0; k<4; ++k)
            sum+= f[k];
        result+= sum/2520.;
        sum= 0.;
        for (int k=0; k<3; ++k)
            sum+= f[EdgeByVert(i,k)+4];
        result+= -sum/630.;
        sum= 0.;
        const int oppF=OppFace(i);
        for (int k=0; k<3; ++k)
            sum+= f[EdgeOfFace(oppF, k)+4];
        result+= -sum/420.;
        return result;
    }
    else  // hat function on edge
    {
        i-= 4;
        // Q = sum c[i]*f[i]
        // Gewichte c[i] = 4/315 	fuer egde i
        //                 1/315	fuer opposite edge
        //                 2/315	fuer uebrige edges
        //                -1/630    	fuer an edge i anliegende verts
        //                -1/420        fuer uebrige zwei verts   
        result=  f[i+4]*4./315.;
        const int opp= OppEdge(i);
        result+= f[opp+4]/315.;
        sum= 0.;
        for(int k=0; k<6; ++k)
            if (k!=i && k!=opp)
                sum+= f[k+4];
        result+= sum*2./315.;
        sum= f[VertOfEdge(i,0)] + f[VertOfEdge(i,1)];
        result+= -sum/630.;
        sum= f[VertOfEdge(OppEdge(i),0)] + f[VertOfEdge(OppEdge(i),1)];
        result+= -sum/420.;

        return result;
    }
}
*/

inline double GetMassP2( int i, int j)
{
// zur Erlaeuterung der zurueckgegebenen Gewichte
// beachte Kommentare zu Funktion QuadP2 obendrueber!

    if (i>j) { int h=i; i=j; j=h; } // swap such that i<=j holds
    if (j<4) // i,j are vertex-indices
    {
        if (i==j)
            return 1./420.;
        else
            return 1./2520.;
    }
    else if (i>=4) // i,j are edge-indices
    {
        if (i==j)
            return 4./315.;
        else
        {
            if (i==OppEdge(j))
                return 1./315.;
            else
                return 2./315.;
        }
    }
    else // i<4, j>=4: Sonderfall...
    {
        if (i==VertOfEdge(j,0) || i==VertOfEdge(j,1))
            return -1./630.;
        else
            return -1./420.;
    }
}

template<class DiscVelSolT>
void LevelsetP2CL::SetupSystem( const DiscVelSolT& vel)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> E(&_E, num_unks, num_unks), 
                               H(&_H, num_unks, num_unks);
    IdxT Numb[10];
    
    std::cerr << "entering SetupSystem: " << num_unks << " levelset unknowns. ";

    // fill value part of matrices
    Quad2CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad2CL<double> u_Grad[10]; // fuer u grad v_i
    SMatrixCL<3,3> T;
    double det, absdet, h_T;

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        h_T= pow( absdet, 1./3.);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and u_loc
        for(int i=0; i<4; ++i)
        {
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
	    u_loc[i]= vel.val( *sit->GetVertex(i));
        }
        for(int i=0; i<6; ++i)
        {
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }
        u_loc[4]= vel.val( *sit, 0.25, 0.25, 0.25);

        for(int i=0; i<10; ++i)
            u_Grad[i]= dot( u_loc, Grad[i]);
        
        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            for(int j=0; j<10; ++j)
            {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
                E( Numb[i], Numb[j])+= GetMassP2(i,j) * absdet
                                     + u_Grad[i].quadP2(j, absdet)*_SD*h_T; 
                
                // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                H( Numb[i], Numb[j])+= u_Grad[j].quadP2(i, absdet)
                                     + Quad2CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * _SD*h_T;
            }
    }
    
    E.Build();
    H.Build();
    std::cerr << _E.num_nonzeros() << " nonzeros in E, "
              << _H.num_nonzeros() << " nonzeros in H! " << std::endl;
}

void LevelsetP2CL::Reparam( Uint steps, double dt)
// Reparametrization of the levelset function Phi
{
    VectorCL Psi= Phi.Data, b;
    MatrixCL L, R, M;

    for (Uint i=0; i<steps; ++i)
    {
        SetupReparamSystem( M, R, Psi, b);
        L.LinComb( 1., M, dt*_theta, R);
        
        b*= dt;
        b+= M*Psi - dt*(1.-_theta) * (R*Psi);
        _gm.Solve( L, Psi, b);
        std::cout << "Reparam: res = " << _gm.GetResid() << ", iter = " << _gm.GetIter() << std::endl;
    }
    
    Phi.Data= Psi;
}

double func_abs( double x) { return std::abs(x); }

void LevelsetP2CL::SetupReparamSystem( MatrixCL& _M, MatrixCL& _R, const VectorCL& Psi, VectorCL& b) const
// M, R, b describe the following terms used for reparametrization:  
// b_i  = ( S(Phi0),           v_i     )
// M_ij = ( v_j,               v_i     )
// R_ij = ( w(Psi) grad v_j,   v_i     )
//      + (|S(Phi0)| grad v_j, grad v_i) * diff
// where v_i, v_j denote the ansatz functions
// and w(Psi) = sign(Phi0) * grad Psi / |grad Psi| the scaled gradient of Psi
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> R(&_R, num_unks, num_unks);
    SparseMatBuilderCL<double> M(&_M, num_unks, num_unks);
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
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        
        for (int i=0; i<4; ++i)
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        for (int i=0; i<6; ++i)
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
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
                M( Numb[i], Numb[j])+= GetMassP2(i,j)*absdet;
                // R_ij = ( w(Psi) grad v_j, v_i + SD * w(Psi) grad v_i )
                R( Numb[i], Numb[j])+= w_Grad[j].quadP2(i, absdet)
                    + _diff*Quad2CL<>(dot( Grad[j]*Sign_Phi.apply( func_abs), Grad[i])).quad( absdet);
            }
        }
    }
    M.Build();
    R.Build();
    std::cerr << _M.num_nonzeros() << " nonzeros in M, " 
              << _R.num_nonzeros() << " nonzeros in R!" << std::endl;
}

void LevelsetP2CL::SetTimeStep( double dt, double theta) 
{ 
    _dt= dt; 
    if (theta >= 0) _theta= theta;
    
    _L.LinComb( 1., _E, _theta*_dt, _H); 
}

void LevelsetP2CL::ComputeRhs( VectorCL& rhs) const
{
    rhs= _E*Phi.Data - _dt*(1-_theta) * (_H*Phi.Data);
}

void LevelsetP2CL::DoStep( const VectorCL& rhs)
{
    _gm.Solve( _L, Phi.Data, rhs);
    std::cout << "res = " << _gm.GetResid() << ", iter = " << _gm.GetIter() <<std::endl;
}

void LevelsetP2CL::DoStep()
{
    VectorCL rhs( Phi.Data.size());
    ComputeRhs( rhs);
    DoStep( rhs);
}

bool LevelsetP2CL::Intersects( const TetraCL& t) const
// Teste, ob Phi in allen DoFs das gleiche Vorzeichen hat
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
// Reparametrisierung durch Fast Marching Method
// Dabei verschiebt sich der 0-Level nur innerhalb der Elemente, die diesen schneiden.
// onlyZeroLvl==true  =>  es wird nur lokal an der Phasengrenze reparametrisiert
{
    FastMarchCL fm( _MG, Phi);
    
    if (OnlyZeroLvl)
    {
        if (Periodic)
            fm.InitZeroPer( _Bnd, ModifyZero);
        else
            fm.InitZero( ModifyZero);
        fm.RestoreSigns();
    }
    else if (Periodic)
        fm.ReparamPer( _Bnd, ModifyZero);
    else
        fm.Reparam( ModifyZero);
}


void LevelsetP2CL::AccumulateBndIntegral( VecDescCL& f) const
{
    VecDescCL SmPhi= Phi;
    if (_curvDiff>0)
        SmoothPhi( SmPhi.Data, _curvDiff);
        
    const Uint  idx_f=     f.RowIdx->GetIdx();
    Point2DCL   AT_i, ab, tmp;
    IdxT        Numb[10];
    
//std::ofstream fil("surf.off");
//fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";

    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    InterfacePatchCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
 
    for (MultiGridCL::TriangTetraIteratorCL it=_MG.GetTriangTetraBegin(), end=_MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, tmp[0], *it);
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
                
                Point3DCL gr, gradv[4];
                // gr= grad Phi(P) + grad Phi(Q) + grad Phi(R)
                for (int node=0; node<4; ++node)
                { // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
                    gradv[node]= Grad[v][node];
                    gr+= BaryPQR[node]*gradv[node];
                }
                    
                if (patch.IsQuadrilateral())
                {
                    Point3DCL grSQR;
                    // grSQR = grad Phi(S) + grad Phi(Q) + grad Phi(R)
                    for (int k=0; k<4; ++k)
                        grSQR+= BarySQR[k]*gradv[k];

                    gr+= patch.GetAreaFrac() * grSQR;
                }
                // nun gilt:
                // gr = [grad Phi(P)+...] + (a+b-1)[grad Phi(S)+...]
                
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

double LevelsetP2CL::GetVolume( double translation) const
{
    const Uint lvl= Phi.GetLevel();
    const_DiscSolCL phi= GetSolution();
    SmoothedJumpCL H( JumpCL( 1, 0), DROPS::H_sm, 1e-4);
    Quad2CL<> Xi;    // 1 fuer phi<0, 0 sonst
    double det, absdet, vol= 0;
    SMatrixCL<3,3> T;
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        absdet= fabs( det);
        
        for (int i=0; i<4; ++i)
            Xi[i]= H( phi.val( *sit->GetVertex(i)) + translation);
        // values in barycenter
        Xi[4]= H( phi.val( *sit, 0.25, 0.25, 0.25) + translation);
        
        vol+= Xi.quad( absdet);
    }
    return vol;
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
        const double d2=1.2*(v1*d0-v0*d1)/(v1-v0);
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
    PCG_SsorCL pcg( _pc, _gm.GetMaxIter(),_gm.GetTol());
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
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> Mb(&M, num_unks, num_unks);
    SparseMatBuilderCL<double> Ab(&A, num_unks, num_unks);

    Quad2CL<Point3DCL> Grad[10], GradRef[10], w_loc;
    SMatrixCL<3,3>     T;
    P2DiscCL::GetGradientsOnRef( GradRef);
    
    IdxT         Numb[10];
    double det, absdet;
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        
        for (int i=0; i<4; ++i)
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        for (int i=0; i<6; ++i)
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            for(int j=0; j<10; ++j)
            {
                // M_ij = ( v_j, v_i),    A_ij = ( grad v_j, grad v_i)
                Mb( Numb[i], Numb[j])+= GetMassP2(i,j)*absdet;
                Ab( Numb[i], Numb[j])+= Quad2CL<>(dot( Grad[j], Grad[i])).quad( absdet);
            }
        }
    }
    Mb.Build();
    Ab.Build();
}

InterfacePatchCL::InterfacePatchCL() 
  : approxZero_(1e-8), RegRef_( GetRefRule( RegRefRuleC)), intersec_(0)
{
    BaryDoF_[0][0]= BaryDoF_[1][1]= BaryDoF_[2][2]= BaryDoF_[3][3]= 1.;
    for (int edge=0; edge<6; ++edge)
        BaryDoF_[edge+4]= 0.5*(BaryDoF_[VertOfEdge(edge,0)] + BaryDoF_[VertOfEdge(edge,1)]);
}

void InterfacePatchCL::Init( TetraCL& t, VecDescCL& ls)
{
    const Uint idx_ls= ls.RowIdx->GetIdx();
    for (int v=0; v<10; ++v)
    { // collect data on all DoF
        Coord_[v]= v<4 ? t.GetVertex(v)->GetCoord() : GetBaryCenter( *t.GetEdge(v-4));
        PhiLoc_[v]= ls.Data[v<4 ? t.GetVertex(v)->Unknowns(idx_ls) : t.GetEdge(v-4)->Unknowns(idx_ls)];
        sign_[v]= std::abs(PhiLoc_[v]) < approxZero_ ? 0 : (PhiLoc_[v]>0 ? 1 : -1);
    }
}

bool InterfacePatchCL::ComputeForChild( Uint ch)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch]);
    num_sign_[0]= num_sign_[1]= num_sign_[2]= 0;
    for (int vert= 0; vert<4; ++vert)
        ++num_sign_[ sign_[data.Vertices[vert]] + 1];

    intersec_= 0;
    if (num_sign_[0]*num_sign_[2]==0 && num_sign_[1]<3) // no change of sign on child
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
            const double lambda= PhiLoc_[v0]/(PhiLoc_[v0]-PhiLoc_[v1]);
            Bary_[intersec_]= (1-lambda)*BaryDoF_[v0] + lambda * BaryDoF_[v1];
            // bary-coords of tetra, not of subtetra!
            PQRS_[intersec_++]= (1-lambda) * Coord_[v0] + lambda * Coord_[v1];
        }
    }
    if (intersec_<3) return false; // Nullstellenmenge vom Mass 0!

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
        // Loese (Q-P)a + (R-P)b = S-P
        SMatrixCL<2,2> M;  
        M(0,0)= A(0,0); M(0,1)= A(0,1);          // 1st row of A
        int row2= 1;
	double detM= A(0,0)*A(1,1) - A(1,0)*A(0,1);
        if (std::abs( detM)<1e-15 ) // upper 2x2 part of A close to singular
        {
            row2= 2;
            detM= A(0,0)*A(2,1) - A(2,0)*A(0,1);
        }
        M(1,0)= A(row2,0); M(1,1)= A(row2,1);
        // now M is nonsingular 2x2 part of A
        tmp[0]= PQRS_[3][0]-PQRS_[0][0]; tmp[1]= PQRS_[3][row2]-PQRS_[0][row2];
        // tmp = S-P
        Solve2x2( detM, M, ab_, tmp);
        //if (ab_[0]<0 || ab_[1]<0) 
        //    std::cerr<<"LevelsetP2CL::AccumulateBndIntegral: a or b negative"<<std::endl;
        // a,b>=0 muss erfuellt sein, da wegen edge+oppEdge==5 die Punkte P und S sich automatisch gegenueber liegen muessten...
    }

    if (EqualToFace()) // interface is shared by two tetras
        sqrtDetATA_/= 2;
        
    return true; // computed patch of child;
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

} // end of namespace DROPS


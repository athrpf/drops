//**************************************************************************
// File:    levelset.tpp                                                   *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include <fstream>

namespace DROPS
{

inline double Sign( double x)
{
    return x<0 ? -1 : x>0 ? 1 : 0;
}

inline double SmoothedSign( double x, double alpha)
{
    return x/sqrt(x*x+alpha*alpha);
}

void LevelsetP2CL::CreateNumbering(Uint level, IdxDescCL* idx)
// used for numbering of the Unknowns depending on the index IdxDesc[idxnum].
// sets up the description of the index idxnum in IdxDesc[idxnum],
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, that IdxDesc[idxnum].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx->TriangLevel = level;
    idx->NumUnknowns = 0;

    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL

    // allocate space for indices; number unknowns in TriangLevel level
    CreateNumbOnVertex( idxnum, idx->NumUnknowns, idx->NumUnknownsVertex,
                        _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level),
                        _dummyBnd);
    CreateNumbOnEdge( idxnum, idx->NumUnknowns, idx->NumUnknownsEdge,
                      _MG.GetTriangEdgeBegin(level), _MG.GetTriangEdgeEnd(level),
                      _dummyBnd );
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
    const Uint lvl= Phi.RowIdx->TriangLevel,
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
    const Uint lvl= Phi.RowIdx->TriangLevel,
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
	    u_loc.val[i]= vel.val( *sit->GetVertex(i));
        }
        for(int i=0; i<6; ++i)
        {
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }
        u_loc.val[4]= vel.val( *sit, 0.25, 0.25, 0.25);

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
    MatrixCL L, R;

    for (Uint i=0; i<steps; ++i)
    {
        SetupReparamSystem( R, Psi, b);
        L.LinComb( 1., _E, dt*_theta, R);
        
        b*= dt;
        b+= _E*Psi - dt*(1.-_theta) * (R*Psi);
        _gm.Solve( L, Psi, b);
        std::cout << "Reparam: res = " << _gm.GetResid() << ", iter = " << _gm.GetIter() << std::endl;
    }
    
    Phi.Data= Psi;
}

double func_abs( double x) { return std::abs(x); }

void LevelsetP2CL::SetupReparamSystem( MatrixCL& _R, const VectorCL& Psi, VectorCL& b)
// R, b describe the following terms used for reparametrization:  
// b_i  = ( S(Phi0),           v_i              + SD * w(Psi) grad v_i )
// R_ij = ( w(Psi) grad v_j,   v_i              + SD * w(Psi) grad v_i )
//      + (|S(Phi0)| grad v_j, grad v_i) * diff
// where v_i, v_j denote the ansatz functions
// and w(Psi) = sign(Phi0) * grad Psi / |grad Psi| the scaled gradient of Psi
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.RowIdx->TriangLevel,
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> R(&_R, num_unks, num_unks);
    b.resize( 0);
    b.resize( num_unks);
    std::cerr << "entering SetupReparamSystem: " << num_unks << " levelset unknowns. ";

    Quad2CL<double>    Sign_Phi;
    Quad2CL<Point3DCL> Grad[10], GradRef[10], w_loc;
    SMatrixCL<3,3>     T;
    P2DiscCL::GetGradientsOnRef( GradRef);
    
    IdxT         Numb[10];
    SVectorCL<3> grad_Psi[4];
    DiscSolCL    phi= GetSolution();
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
            Sign_Phi.val[i]= SmoothedSign( phi.val( *sit->GetVertex(i)), alpha);
            grad_Psi[i]= Point3DCL();  // init with zero
            for (int l=0; l<10; ++l)
                grad_Psi[i]+= Psi[ Numb[l]] * Grad[l].val[i];
	    w_loc.val[i]= (Sign_Phi.val[i]/grad_Psi[i].norm() )*grad_Psi[i];
        }
        // values in barycenter
        Point3DCL gr= 0.25*(grad_Psi[0]+grad_Psi[1]+grad_Psi[2]+grad_Psi[3]);
        Sign_Phi.val[4]= SmoothedSign( phi.val( *sit, 0.25, 0.25, 0.25), alpha);
        w_loc.val[4]=    Sign_Phi.val[4]*gr/gr.norm();

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
        {
            // b_i  = ( S(Phi0),         v_i + SD * w(Psi) grad v_i )
            b[ Numb[i]]+= Sign_Phi.quadP2( i, absdet);
//            b[ Numb[i]]+= _SD*h_T*QuadVelGrad(w_loc,Grad[i], Sign_Phi)*absdet; 
            for(int j=0; j<10; ++j)
            {
                // R_ij = ( w(Psi) grad v_j, v_i + SD * w(Psi) grad v_i )
                R( Numb[i], Numb[j])+= Quad2CL<>(dot( w_loc, Grad[j])).quadP2(i, absdet)
                    + _diff*Quad2CL<>(dot( Grad[j]*Sign_Phi.apply( func_abs), Grad[i])).quad( absdet);
//                R( Numb[i], Numb[j])+= _SD*h_T*QuadVelGrad(w_loc,Grad[j],Grad[i])*absdet;
            }
        }
    }
    
    R.Build();
    std::cerr << _R.num_nonzeros() << " nonzeros in R!" << std::endl;
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

/*
void LevelsetP2CL::Reparam2()
// Reparametrisierung, die den 0-Levelset fest laesst
// Idee: Phi auf T durch festes g=|grad Phi|>0 teilen -> 0-Levelset bleibt fest
{
    const size_t num_unks= Phi.RowIdx->NumUnknowns;
    const Uint   idx=      Phi.RowIdx->GetIdx();
    VectorBaseCL<bool> nearby( false, num_unks);  // is DoF near the levelset?
    VectorBaseCL<int>  num( 0, num_unks);         // how many tetras have this DoF?
    VectorCL           div( num_unks);            // sum of |grad Phi| on this tetras
    int                num_ls= 0;                 // how many tetras intersect the levelset?
    double             div_ls= 0;                 // sum of |grad Phi| on this tetras
    double      PhiLoc[10], det;
    int         sign[10];
    int         num_sign[3]; // - 0 + 
    IdxT        Numb[10];
    Point3DCL   grad;
    
    SMatrixCL<3,3> T;
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    for (MultiGridCL::TriangTetraIteratorCL it=_MG.GetTriangTetraBegin(), end=_MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            if (v<4) 
            {
                PhiLoc[v]= Phi.Data[it->GetVertex(v)->Unknowns(idx)];
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            }
            else
            {
                PhiLoc[v]= Phi.Data[it->GetEdge(v-4)->Unknowns(idx)];
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);
            }
            sign[v]= std::abs(PhiLoc[v])<1e-8 ? 0 : (PhiLoc[v]>0 ? 1 : -1);
        }
            
        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];
            const bool intersec= !(num_sign[0]*num_sign[2]==0 && num_sign[1]!=3); // (...) = no change of sign on child
            const double Phi0= PhiLoc[data.Vertices[0]];
            for (int i=0; i<3; ++i)
                grad[i]= PhiLoc[data.Vertices[i+1]] - Phi0;
            const double g= (T*grad).norm()*2; // doppelt so gross, da hier Kind des Tetra betrachtet wird

            if (intersec)
            {
                div_ls+= g;    ++num_ls;
                for (int vert= 0; vert<4; ++vert)
                    nearby[Numb[data.Vertices[vert]]]= true;
            }
            else
                for (int vert= 0; vert<4; ++vert)
                {
                    const size_t nv= Numb[data.Vertices[vert]];
                    div[nv]+= g;    ++num[nv];
                }
        }
    }
    
    const double c= num_ls/div_ls;
    std::cerr << "ReparamSaveIF: " << num_ls << " child tetras intersecting the levelset\n"
              << "\twith average norm of gradient = " << 1./c << std::endl;
    
    for (size_t i=0; i<num_unks; ++i)
        if (nearby[i])
            Phi.Data[i]*= c;
        else if (div[i]>1e-4)
        {
const double d= div[i]/num[i];
// Ausgabe nur bei rel. Abweichung von >50%
if (fabs(d-1./c)*c>0.5) std::cerr << d << "\t";
            Phi.Data[i]*= num[i]/div[i];
        }
std::cerr << std::endl;        
}
*/

inline void Solve2x2( const SMatrixCL<2,2>& A, SVectorCL<2>& x, const SVectorCL<2>& b)
{
    const double det= A(0,0)*A(1,1) - A(0,1)*A(1,0);
    x[0]= (A(1,1)*b[0]-A(0,1)*b[1])/det;
    x[1]= (A(0,0)*b[1]-A(1,0)*b[0])/det;
}

void LevelsetP2CL::AccumulateBndIntegral( VecDescCL& f) const
{
    BaryCoordCL BaryDoF[10];
    Point3DCL   Coord[10];
    BaryDoF[0][0]= BaryDoF[1][1]= BaryDoF[2][2]= BaryDoF[3][3]= 1.;
    for (int edge=0; edge<6; ++edge)
        BaryDoF[edge+4]= 0.5*(BaryDoF[VertOfEdge(edge,0)] + BaryDoF[VertOfEdge(edge,1)]);

    const Uint  idx_phi= Phi.RowIdx->GetIdx(),
                idx_f=     f.RowIdx->GetIdx();
    double      PhiLoc[10];
    int         sign[10];
    int         num_sign[3]; // - 0 + 
    BaryCoordCL Bary[4];
    Point3DCL   PQRS[4], B[3];
    Point2DCL   AT_i, ab, tmp;
    IdxT        Numb[10];
/*    
std::ofstream fil("surf.off");
fil << "appearance {\n-concave\nshading smooth\n}\nLIST\n{\n";
*/
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    RefRuleCL RegRef= GetRefRule( RegRefRuleC);

    P2DiscCL::GetGradientsOnRef( GradRef);
 
    for (MultiGridCL::TriangTetraIteratorCL it=_MG.GetTriangTetraBegin(), end=_MG.GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, tmp[0], *it);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder

        for (int v=0; v<10; ++v)
        { // collect data on all DoF
            Coord[v]= v<4 ? it->GetVertex(v)->GetCoord() : GetBaryCenter( *it->GetEdge(v-4));
            if (v<4) 
            {
                PhiLoc[v]= Phi.Data[it->GetVertex(v)->Unknowns(idx_phi)];
                Numb[v]= it->GetVertex(v)->Unknowns.Exist(idx_f) ?
                            it->GetVertex(v)->Unknowns(idx_f) : NoIdx;
            }
            else
            {
                PhiLoc[v]= Phi.Data[it->GetEdge(v-4)->Unknowns(idx_phi)];
                Numb[v]= it->GetEdge(v-4)->Unknowns.Exist(idx_f) ?
                             it->GetEdge(v-4)->Unknowns(idx_f) : NoIdx;
            }
            sign[v]= std::abs(PhiLoc[v])<1e-8 ? 0 : (PhiLoc[v]>0 ? 1 : -1);
        }
            
        for (int ch=0; ch<8; ++ch)
        {
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for (int vert= 0; vert<4; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];
            if (num_sign[0]*num_sign[2]==0 && num_sign[1]!=3) // no change of sign on child
                continue;
            if (num_sign[1]==4)
            { 
                std::cerr << "WARNING: LevelsetP2CL::AccumulateBndIntegral: found 3-dim. zero level set, grid is too coarse!" << std::endl; 
                continue; 
            }

            int intersec= 0;
            // erst werden die Nullknoten in PQRS gespeichert...
            for (int vert= 0; vert<4; ++vert)
            {
                const int v= data.Vertices[vert];
                if (sign[v]==0)
                {
                    Bary[intersec]= BaryDoF[v];
                    PQRS[intersec++]= Coord[v];
                }
            }
            // ...dann die echten Schnittpunkte auf den Kanten mit Vorzeichenwechsel
            for (int edge= 0; edge<6; ++edge)
            {
                const int v0= data.Vertices[ VertOfEdge( edge, 0)],
                          v1= data.Vertices[ VertOfEdge( edge, 1)];
                if (sign[v0]*sign[v1]<0) // different sign -> 0-level intersects this edge
                {
                    const double lambda= PhiLoc[v0]/(PhiLoc[v0]-PhiLoc[v1]);
                    Bary[intersec]= (1-lambda)*BaryDoF[v0] + lambda * BaryDoF[v1];
                    // bary-coords of tetra, not of subtetra!
                    PQRS[intersec++]= (1-lambda) * Coord[v0] + lambda * Coord[v1];
                }
            }
/*
fil << "geom {OFF " << intersec << " 1 0\n";
for (int i=0; i<intersec; ++i)
{
    for (int j=0; j<3; ++j)
        fil << PQRS[i][j] << ' ';
    fil << '\n';
}
if (intersec==4)
    fil << "4 0 1 3 2";
else
    fil << "3 0 1 2";
fil << "\n}\n";
*/

            if (intersec<3) continue; // Nullstellenmenge vom Mass 0!

            SMatrixCL<3,2> A;    // A = [ Q-P | R-P ]
            A(0,0)= PQRS[1][0]-PQRS[0][0];    A(0,1)= PQRS[2][0]-PQRS[0][0];
            A(1,0)= PQRS[1][1]-PQRS[0][1];    A(1,1)= PQRS[2][1]-PQRS[0][1];
            A(2,0)= PQRS[1][2]-PQRS[0][2];    A(2,1)= PQRS[2][2]-PQRS[0][2];
            SMatrixCL<2,2> ATA; 
            ATA(0,0)=           A(0,0)*A(0,0)+A(1,0)*A(1,0)+A(2,0)*A(2,0);
            ATA(0,1)= ATA(1,0)= A(0,0)*A(0,1)+A(1,0)*A(1,1)+A(2,0)*A(2,1);
            ATA(1,1)=           A(0,1)*A(0,1)+A(1,1)*A(1,1)+A(2,1)*A(2,1);
            double sqrtDetATA= std::sqrt( ATA(0,0)*ATA(1,1) - ATA(1,0)*ATA(1,0) );
            BaryCoordCL BaryPQR, BarySQR;
            for (int i=0; i<3; ++i)
            {
                // addiere baryzentrische Koordinaten von P,Q,R bzw. S,Q,R
                BaryPQR+= Bary[i];
                BarySQR+= Bary[i+1];
            
                // berechne B = A * (ATA)^-1 * AT
                AT_i[0]= A(i,0); AT_i[1]= A(i,1);
                Solve2x2(ATA,tmp,AT_i);
                B[i]= A*tmp;
            }

            if (intersec==4) // 4 intersections --> a+b != 1
            { // berechne a, b
                // Loese (Q-P)a + (R-P)b = S-P
                SMatrixCL<2,2> M;  
                M(0,0)= A(0,0); M(0,1)= A(0,1); M(1,0)=A(1,0); M(1,1)= A(1,1);
                // M is upper 2x2 part of A
                tmp[0]= PQRS[3][0]-PQRS[0][0]; tmp[1]= PQRS[3][1]-PQRS[0][1];
                // tmp = S-P
                Solve2x2( M, ab, tmp);

                //if (ab[0]<0 || ab[1]<0) 
                //    std::cerr<<"LevelsetP2CL::AccumulateBndIntegral: a or b negative"<<std::endl;
                // a,b>=0 muss erfuellt sein, da wegen edge+oppEdge==5 die Punkte P und S sich automatisch gegenueber liegen muessten...
            }

            sqrtDetATA/= 6;

            
            for (int v=0; v<10; ++v)
            {
                if (Numb[v]==NoIdx) continue;
                
                Point3DCL gr, gradv[4];
                // gr= grad Phi(P) + grad Phi(Q) + grad Phi(R)
                for (int node=0; node<4; ++node)
                { // gradv = Werte von grad Hutfunktion fuer DoF v in den vier vertices
                    gradv[node]= Grad[v].val[node];
                }
                    
                for (int k=0; k<4; ++k)
                    gr+= BaryPQR[k]*gradv[k];

                if (intersec==4)
                {
                    Point3DCL grSQR;
                    // grSQR = grad Phi(S) + grad Phi(Q) + grad Phi(R)
                    for (int k=0; k<4; ++k)
                        grSQR+= BarySQR[k]*gradv[k];

                    gr+= (ab[0] + ab[1] - 1) * grSQR;
                }
                // nun gilt:
                // gr = [grad Phi(P)+...] + (a+b-1)[grad Phi(S)+...]
                
                for (int i=0; i<3; ++i)
                {
                    const double val= inner_prod( gr, B[i]);
                    f.Data[Numb[v]+i]-= sigma * sqrtDetATA * val;
                }
            }
        } // Ende der for-Schleife ueber die Kinder
    }
//fil << "}\n";    
}


} // end of namespace DROPS


/// \file
/// \brief Classes that constitute a 2-phase-transport-problem.
/// change Quad5CL (15 points) to Quad3CL (5 points) in Setup matrices routines

#include "transport/transport2phase.h"
#include "transport/localsetups.cpp"

namespace DROPS
{
//=====================================
//TransportP1CL
//=====================================
void TransportP1CL::Init (instat_scalar_fun_ptr cneg, instat_scalar_fun_ptr cpos, double t)
{
    ct.t = t;
    const Uint lvl= ct.GetLevel(),
               idx= ct.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( idx)) {
            if (lset_.Phi.Data[it->Unknowns( lset_.Phi.RowIdx->GetIdx())] <= 0.)
                ct.Data[it->Unknowns( idx)]= H_*cneg( it->GetCoord(), t);
            else
                ct.Data[it->Unknowns( idx)]= cpos( it->GetCoord(), t);
        }
    }
}

void TransportP1CL::SetTimeStep (double dt, double theta)
{
    dt_= dt;
    if (theta >= 0. && theta <= 1.) theta_= theta;
}

void TransportP1CL::SetupInstatSystem (MatrixCL& matA, VecDescCL* cplA,
                        MatrixCL& matM, VecDescCL* cplM, MatrixCL& matC, VecDescCL* cplC, VecDescCL* cplb, VecDescCL* oldcplM,
                        IdxDescCL& RowIdx, const double time, const double time_old) const
{
    if (cplM != 0)
    {
        cplM->Data= 0.;
        oldcplM->Data= 0.;
        cplb->Data= 0.;
        cplA->Data= 0.;
        cplC->Data= 0.;
    }
    matM.clear();
    matA.clear();
    matC.clear();
    const IdxT num_unks=  RowIdx.NumUnknowns();
    MatrixBuilderCL A(&matA, num_unks,  num_unks),//diffusion
                    M(&matM, num_unks,  num_unks),//mass matrix
                    C(&matC, num_unks,  num_unks);// convection

    const Uint lvl= RowIdx.TriangLevel();
    LocalNumbP1CL n;
    double coupA[4][4]={0.}, coupM[4][4]={0.}, coupC[4][4]={0.}, locf[4]={0.};

    // The 16 products of the P1-shape-functions
    LocalP1CL<> p1[4];
    Quad3CL<> q3_p[4];
    Quad5CL<> q5_p[4];
    for(int i= 0; i < 4; ++i) {
        p1[i][i]=1.;
        q3_p[i].assign(p1[i]);
        q5_p[i].assign(p1[i]);
     }
    LocalP2CL<> pipj[4][4];
    SetupPiPj(pipj);
    const P2EvalCL<SVectorCL<3>, const VelBndDataT, VecDescCL> u( v_, &Bnd_v_, &MG_, time);

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {
        n.assign( *sit, RowIdx, Bnd_);      
        double det, absdet;
        SMatrixCL<3,4> G;
        P1DiscCL::GetGradients( G, det, *sit);
        const SMatrixCL<4,4> GTG( GramMatrix( G));
        absdet= std::fabs( det);  
        Quad5CL<> rhs( *sit, f_, time);
        Quad3CL<Point3DCL> q3_u( *sit, u, time);
        LocalP2CL<Point3DCL> lp2_u( *sit, u, time);
        InterfacePatchCL cut;
        cut.Init( *sit, lset_.Phi);
        SetupLocalRhs( locf, rhs, q5_p, absdet);        
        bool nocut=!cut.Intersects();
        if (nocut)
        {
            bool pPart= (cut.GetSign( 0) == 1);
            SetupLocalOnePhaseSystem ( coupM, coupA, coupC, q3_u, q3_p, G, GTG, absdet, D_, H_, pPart);
        }
        else
            SetupLocalOneInterfaceSystem( cut, coupM, coupA, coupC, absdet, D_, H_, lp2_u, pipj, p1, G, GTG);
         
        // write values into matrix
        for(int i= 0; i < 4; ++i)
            if (n.WithUnknowns( i)){
                for(int j= 0; j < 4; ++j)
                    if (n.WithUnknowns( j)) {
                        M( n.num[i], n.num[j])+= coupM[j][i];
                        A( n.num[i], n.num[j])+= coupA[j][i];
                        C( n.num[i], n.num[j])+= coupC[j][i];
                    }
                    else if (cplM != 0) {
                        const double val= Bnd_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time),
                                     val_old= Bnd_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time_old);
                        cplM->Data[n.num[i]]-= coupM[j][i]*val;
                        oldcplM->Data[n.num[i]]-= coupM[j][i]*val_old;
                        cplA->Data[n.num[i]]-= coupA[j][i]*val;
                        cplC->Data[n.num[i]]-= coupC[j][i]*val;
                    }
                if (cplb!=0) cplb->Data[n.num[i]]+= locf[i];
            }        
    }
    A.Build();
    M.Build();
    C.Build();
}

void TransportP1CL::SetupInstatSystem (MLMatDescCL& matA, VecDescCL& cplA,
    MLMatDescCL& matM, VecDescCL& cplM, MLMatDescCL& matC, VecDescCL& cplC, VecDescCL& cplb, VecDescCL& oldcplM, const double time, const double time_old) const
{
    MLMatrixCL::iterator itA = --matA.Data.end();
    MLMatrixCL::iterator itM = --matM.Data.end();
    MLMatrixCL::iterator itC = --matC.Data.end();
    MLIdxDescCL::iterator it = --matA.RowIdx->end();
    SetupInstatSystem (*itA, &cplA, *itM, &cplM, *itC, &cplC, &cplb, &oldcplM, *it, time, time_old);
    for (size_t num= 1; num < matA.Data.size(); ++num, --itA, --itM, --itC, --it)
        SetupInstatSystem (*itA, 0, *itM, 0, *itC, 0, 0, 0, *it, time, time_old);
}


void TransportP1CL::Update()
{
    MLIdxDescCL* cidx= &idx;

    c.SetIdx( cidx);
    ct2c();

    cplM.SetIdx( cidx);
    cplA.SetIdx( cidx);
    cplC.SetIdx( cidx);
    b.SetIdx( cidx);

    oldcplM.SetIdx( cidx);
    oldcplA.SetIdx( cidx);
    oldcplC.SetIdx( cidx);
    oldb.SetIdx( cidx);

    M.Data.clear();
    M.SetIdx( cidx, cidx);
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    C.Data.clear();
    C.SetIdx( cidx, cidx);

    SetupInstatSystem( A, oldcplA, M, cplM, C, oldcplC, oldb, oldcplM, t_, t_-dt_);
}

void TransportP1CL::InitStep (VectorCL& rhs)
{
    VectorCL rhs1(ct.Data);
    rhs.resize(ct.Data.size());
    if (theta_ != 1.0) {
        VectorCL tmp( dt_*(1. - theta_)*(oldb.Data + oldcplA.Data - A.Data*ct.Data + oldcplC.Data - C.Data*ct.Data));
        VectorCL rhs2;
        rhs2.resize(ct.Data.size());
        gm_.Solve( M.Data, rhs2, tmp);
        std::cerr << "Inverting M_old: res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
        rhs1+= rhs2;
    }
    SetupInstatSystem( A, cplA, M, cplM, C, cplC, b, oldcplM, t_, t_-dt_);
    // Todo: Add M_new(c_{dirichlet,old}), if there are time-dependant Dirichlet-BC.
    rhs= M.Data*rhs1 + cplM.Data - oldcplM.Data  /*Todo: add together with the above cplM.Data + */ + dt_*theta_*(b.Data + cplA.Data + cplC.Data);

    MLMatrixCL L;
    L.LinComb( theta_, A.Data, theta_, C.Data);
    L_.clear();
    L_.LinComb( 1., M.Data, dt_, L);
}

void TransportP1CL::DoStep (const VectorCL& rhs)
{
    std::cerr << "Before solve: res = " << norm( L_*ct.Data - rhs) << std::endl;
    gm_.Solve( L_, ct.Data, rhs);
    L_.clear();
    std::cerr << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void TransportP1CL::CommitStep ()
{
//     std::swap( cplM, oldcplM);
    std::swap( cplA, oldcplA);
    std::swap( cplC, oldcplC);
    std::swap( b, oldb);
    ct2c();
}

void TransportP1CL::DoStep (double new_t)
{
    VectorCL rhs( c.Data.size());
    t_= new_t;
    InitStep( rhs);
    DoStep( rhs);
    CommitStep();
}

void TransportP1CL::c2ct()
// For lset <= 0, multiply with H
{
    Uint cidx= c.RowIdx->GetIdx(),
         phiidx= lset_.Phi.RowIdx->GetIdx();
    ct.Data= c.Data;
    DROPS_FOR_TRIANG_VERTEX( MG_, c.RowIdx->TriangLevel(), sit)
        if (lset_.Phi.Data[sit->Unknowns( phiidx)] <= 0. && sit->Unknowns.Exist( cidx))
            ct.Data[sit->Unknowns( cidx)]*= H_;
}

void TransportP1CL::ct2c()
// For lset <= 0, divide by H
{
    Uint ctidx= ct.RowIdx->GetIdx(),
         phiidx= lset_.Phi.RowIdx->GetIdx();
    c.Data= ct.Data;
    DROPS_FOR_TRIANG_VERTEX( MG_, ct.RowIdx->TriangLevel(), sit)
        if (lset_.Phi.Data[sit->Unknowns( phiidx)] <= 0. && sit->Unknowns.Exist( ctidx))
            c.Data[sit->Unknowns( ctidx)]/= H_;
}

//*****************************************************************************
//                               TransportRepairCL
//*****************************************************************************
void
TransportRepairCL::post_refine ()
{
    VecDescCL loc_ct;
    IdxDescCL loc_cidx( P1_FE);
    VecDescCL& ct= c_.ct;
    match_fun match= mg_.GetBnd().GetMatchFun();

    loc_cidx.CreateNumbering( mg_.GetLastLevel(), mg_, c_.GetBndData(), match);
    loc_ct.SetIdx( &loc_cidx);
    RepairAfterRefineP1( c_.GetSolution( ct), loc_ct);

    ct.Clear();
    ct.RowIdx->DeleteNumbering( mg_);
    c_.idx.GetFinest().swap( loc_cidx);
    ct.SetIdx( &c_.idx);
    ct.Data= loc_ct.Data;
}

inline void
  VelTranspRepairCL_P1::post_refine ()
{
    VelVecDescCL loc_v;
    IdxDescCL    loc_vidx(vecP2_FE);
    VelVecDescCL& v= v_;
    Uint LastLevel= mg_.GetLastLevel();
    match_fun match= mg_.GetBnd().GetMatchFun();
    loc_vidx.CreateNumbering( mg_.GetLastLevel(), mg_, Bnd_v_, match);
    if (LastLevel != v.RowIdx->TriangLevel()) {
        std::cout << "LastLevel: " << LastLevel
                  << " old v->TriangLevel: " << v.RowIdx->TriangLevel() << std::endl;
        throw DROPSErrCL( "VelocityRepairCL::post_refine: Sorry, not yet implemented.");
    }
    loc_v.SetIdx( &loc_vidx);
    RepairAfterRefineP2( const_DiscVelSolCL( &v, &Bnd_v_, &mg_, c_.t_), loc_v);
    v.Clear();
    (*v.RowIdx).DeleteNumbering( mg_);

    vidx_.swap( loc_vidx);
    v.SetIdx( &vidx_);
    v.Data= loc_v.Data;
}

} // end of namespace DROPS

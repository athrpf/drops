/// \file
/// \brief Classes that constitute a 2-phase-transport-problem.

#include "poisson/transport2phase.h"

namespace DROPS
{

void TransportP1CL::Init (instat_scalar_fun_ptr cneg, instat_scalar_fun_ptr cpos)
{
    const Uint lvl= ct.GetLevel(),
               idx= ct.RowIdx->GetIdx();

    DROPS_FOR_TRIANG_VERTEX( MG_, lvl, it) {
        if (it->Unknowns.Exist( idx))
            if (lset_.Phi.Data[it->Unknowns( lset_.Phi.RowIdx->GetIdx())] <= 0.)
                ct.Data[it->Unknowns( idx)]= H_*cneg( it->GetCoord(), 0.);
            else
                ct.Data[it->Unknowns( idx)]= cpos( it->GetCoord(), 0.);
    }
}

void TransportP1CL::SetTimeStep (double dt, double theta)
{
    dt_= dt;
    if (theta >= 0. && theta <= 1.) theta_= theta;
}


void TransportP1CL::SetupLocalSystem (const TetraCL& t,
    double coupM[4][4], double coupA[4][4], double coupC[4][4], const double time, 
    const LocalP2CL<> p2[4], const LocalP2CL<> pipj[4][4], const Quad5CL<> p[4]) const
{
    static LocalP2CL<> one( 1.);

    P2EvalCL<SVectorCL<3>, const VelBndDataT, VecDescCL> u( v_, &Bnd_v_, &MG_, time);

    double det;
    SMatrixCL<3,4> G;
    P1DiscCL::GetGradients( G, det, t);
    const double absdet= std::fabs( det);
    const SMatrixCL<4,4> GTG( GramMatrix( G));

    InterfacePatchCL cut;
    cut.Init( t, lset_.Phi);
    if (!cut.Intersects()) {
        const Quad5CL<Point3DCL> u_loc( t, u, t_);
        const double coeff_d = cut.GetSign( 0) == 1 ? D_[0] : D_[1];
        const double coeff_h = cut.GetSign( 0) == 1 ? 1. : 1./H_;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < i; ++j) {
                coupM[j][i]= coeff_h*P1DiscCL::GetMass( i, j)*absdet;
                coupM[i][j]= coupM[j][i];
                coupA[j][i]= coeff_h*coeff_d*GTG( i, j)/6.0*absdet;
                coupA[i][j]= coupA[j][i];
                coupC[j][i]= coeff_h*Quad5CL<>( dot( u_loc, Quad5CL<Point3DCL>( G.col( j)))*p[i]).quad( absdet);
                coupC[i][j]= coeff_h*Quad5CL<>( dot( u_loc, Quad5CL<Point3DCL>( G.col( i)))*p[j]).quad( absdet);
            }
            coupM[i][i]= coeff_h*P1DiscCL::GetMass( i, i)*absdet;
            coupA[i][i]= coeff_h*coeff_d*GTG( i, i)/6.0*absdet;
            coupC[i][i]= coeff_h*Quad5CL<>( dot( u_loc, Quad5CL<Point3DCL>( G.col( i)))*p[i]).quad( absdet);
        } 
    }
    else {
        const LocalP2CL<Point3DCL> u_loc( t, u, t_);
        LocalP2CL<> convection_interface;
        double iAp, iAn, iMp, iMn, iCp, iCn, integralp, integraln;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < 4; ++j) {
                integralp= integraln= 0.;
                iAp= iAn = iMp = iMn = iCp = iCn = 0.;
                convection_interface= dot( u_loc, LocalP2CL<Point3DCL>(G.col( j)))*p2[i];
                for (int ch= 0; ch < 8; ++ch) {
                    cut.ComputeCutForChild( ch);
                    cut.quadBothParts( integralp, integraln, pipj[i][j], absdet);
                    iMp+= integralp; iMn+= integraln;
                    cut.quadBothParts( integralp, integraln, one, absdet);
                    iAp+= integralp; iAn+= integraln;
                    cut.quadBothParts( integralp, integraln, convection_interface, absdet);
                    iCp+= integralp; iCn+= integraln;
                }
                coupM[j][i]= iMp + iMn/H_;
                coupA[j][i]= (iAp*D_[0] + iAn*D_[1]/H_)* GTG( j, i);
                coupC[j][i]= iCp + iCn/H_;
            }
        }
    }
}


void TransportP1CL::SetupInstatSystem (MatDescCL& matA, VecDescCL& cplA,
    MatDescCL& matM, VecDescCL& cplM, MatDescCL& matC, VecDescCL& cplC, const double time) const
{
    cplM.Data= 0.;
    cplA.Data= 0.;
    cplC.Data= 0.;
    matM.Data.clear();
    matA.Data.clear();
    matC.Data.clear();
    const IdxT num_unks=  matA.RowIdx->NumUnknowns;
    MatrixBuilderCL A(&matA.Data, num_unks,  num_unks),//diffusion
                    M(&matM.Data, num_unks,  num_unks),//mass matrix
                    C(&matC.Data, num_unks,  num_unks);// convection

    const Uint lvl= matA.GetRowLevel();
    LocalNumbP1CL n;
    double coupA[4][4], coupM[4][4], coupC[4][4];

    // The 16 products of the P1-shape-functions
    LocalP2CL<> p2[4], pipj[4][4], convection_interface;
    Quad5CL<> p[4];
    for(int i= 0; i < 4; ++i) {
        LocalP1CL<> p1;
        p1[i]= 1.;
        p[i].assign( p1);
        p2[i].assign( p1);
    }

    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
            pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
        }
        pipj[i][i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
                pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }

    DROPS_FOR_TRIANG_TETRA( MG_, lvl, sit) {
        SetupLocalSystem ( *sit, coupM, coupA, coupC, time, p2, pipj, p);
        n.assign( *sit, *matM.RowIdx, Bnd_);
        // write values into matrix
        for(int i= 0; i < 4; ++i)
            if (n.WithUnknowns( i))
                for(int j= 0; j < 4; ++j)
                    if (n.WithUnknowns( j)) {
                        M( n.num[i], n.num[j])+= coupM[j][i];
                        A( n.num[i], n.num[j])+= coupA[j][i];
                        C( n.num[i], n.num[j])+= coupC[j][i];
                    }
                    else {
                        const double val= Bnd_.GetBndFun( n.bndnum[j])( sit->GetVertex( j)->GetCoord(), time);
                        cplM.Data[n.num[i]]-= coupM[j][i]*val;
                        cplA.Data[n.num[i]]-= coupA[j][i]*val;
                        cplC.Data[n.num[i]]-= coupC[j][i]*val;
                    }
    }
    A.Build();
    M.Build();
    C.Build();
}

void TransportP1CL::Update()
{
    IdxDescCL* cidx= ct.RowIdx;

    c.SetIdx( cidx);
    ct2c();

    cplM.SetIdx( cidx);
    cplA.SetIdx( cidx);
    cplC.SetIdx( cidx);

    oldcplM.SetIdx( cidx);
    oldcplA.SetIdx( cidx);
    oldcplC.SetIdx( cidx);

    M.Data.clear();
    M.SetIdx( cidx, cidx);
    A.Data.clear();
    A.SetIdx( cidx, cidx);
    C.Data.clear();
    C.SetIdx( cidx, cidx);

    SetupInstatSystem( A, oldcplA, M, oldcplM, C, oldcplC, t_);
}

void TransportP1CL::InitStep (VectorCL& rhs)
{
    rhs= - (1. - theta_)*(C.Data*ct.Data - oldcplC.Data);
    SetupInstatSystem( A, cplA, M, cplM, C, cplC, t_);
    rhs+= (1./dt_)*(M.Data*ct.Data - oldcplM.Data + cplM.Data) + theta_*(cplA.Data + cplC.Data)
        - (1. - theta_)*(A.Data*ct.Data - oldcplA.Data);

    MatrixCL L;
    L.LinComb( theta_, A.Data, theta_, C.Data);
    L_.clear();
    L_.LinComb( 1./dt_, M.Data, 1., L);
}

void TransportP1CL::DoStep (const VectorCL& rhs)
{
    std::cerr << "Before solve: res = " << norm( L_*ct.Data - rhs) << std::endl;
    // std::cout << "ct:\n" << ct.Data << "\ncplC:\n" << cplC.Data << "\nC:\n" << C.Data << std::endl;
    // std::cout << "\ncplM:\n" << cplM.Data << "\nM:\n" << M.Data << std::endl;
    // std::cout << "\ncplA:\n" << cplA.Data << "\nA:\n" << A.Data << std::endl;
    gm_.Solve( L_, ct.Data, rhs);
    std::cerr << "res = " << gm_.GetResid() << ", iter = " << gm_.GetIter() << std::endl;
}

void TransportP1CL::CommitStep ()
{
    std::swap( cplM, oldcplM);
    std::swap( cplA, oldcplA);
    std::swap( cplC, oldcplC);
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
    DROPS_FOR_TRIANG_VERTEX( MG_, c.RowIdx->TriangLevel, sit)
        if (lset_.Phi.Data[sit->Unknowns( phiidx)] <= 0. && sit->Unknowns.Exist( cidx))
            ct.Data[sit->Unknowns( cidx)]*= H_;
}

void TransportP1CL::ct2c()
// For lset <= 0, divide by H
{
    Uint ctidx= ct.RowIdx->GetIdx(),
         phiidx= lset_.Phi.RowIdx->GetIdx();
    c.Data= ct.Data;
    DROPS_FOR_TRIANG_VERTEX( MG_, ct.RowIdx->TriangLevel, sit)
        if (lset_.Phi.Data[sit->Unknowns( phiidx)] <= 0. && sit->Unknowns.Exist( ctidx))
            c.Data[sit->Unknowns( ctidx)]/= H_;
}


} // end of namespace DROPS

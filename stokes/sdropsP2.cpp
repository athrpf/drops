#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include <fstream>


inline DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p)
{
    DROPS::SVectorCL<3> ret;
    ret[0]=    sin(p[0])*sin(p[1])*sin(p[2]);
    ret[1]=  - cos(p[0])*cos(p[1])*sin(p[2]);
    ret[2]= 2.*cos(p[0])*sin(p[1])*cos(p[2]);
    return ret/3.;
}

inline double LsgPr(const DROPS::Point3DCL& p)
{
    return -cos(p[0])*sin(p[1])*sin(p[2]);
//     return 1.;
}

// boundary value functions (in 2D-bnd-coords)
//DROPS::StokesVelBndDataCL::bnd_type bnd_val_e2e3(const Point2DCL&);
//DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e3(const Point2DCL&);
//DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e2(const Point2DCL&);


// q*u - nu*laplace u - Dp = f
//                   div u = 0
class StokesCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p)
        { DROPS::SVectorCL<3> ret(0.0); ret[2]= 3.*cos(p[0])*sin(p[1])*cos(p[2]); return ret; }
/*    {
        SVectorCL<3> ret;
        ret[0]=    sin(p[0])*sin(p[1])*sin(p[2]);
        ret[1]=  - cos(p[0])*cos(p[1])*sin(p[2]);
        ret[2]= 2.*cos(p[0])*sin(p[1])*cos(p[2]);
        return ret;
    }
*/    
    const double nu;
    
    StokesCoeffCL() : nu(1.0) {}
};

namespace DROPS // for Strategy
{

template<class MGB, class Coeff>
void Strategy(StokesP2P1CL<MGB,Coeff>& Stokes, double omega, double inner_iter_tol, Uint maxStep)
// flow control
{
    typedef typename StokesP2P1CL<MGB,Coeff>::VelVecDescCL VelVecDescCL;
    
    MultiGridCL& MG= Stokes.GetMG();

    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= &Stokes.vel_idx;
    IdxDescCL* pidx1= &Stokes.pr_idx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx2= &loc_pidx;
//    IdxDescCL* err_idx= &_err_idx;
    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &Stokes.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &Stokes.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;
//    VecDescCL* err= &_err;
    MatDescCL* A= &Stokes.A;
    MatDescCL* B= &Stokes.B;
    Uint step= 0;
//    bool new_marks;
//    double akt_glob_err;

    vidx1->Set(0, 3, 3, 0, 0);
    vidx2->Set(1, 3, 3, 0, 0);
    pidx1->Set(2, 1, 0, 0, 0);
    pidx2->Set(3, 1, 0, 0, 0);
    
    TimerCL time;
//    err_idx->Set(5, 0, 0, 0, 1);
    do
    {
//        akt_glob_err= glob_err;
        MarkAll(MG); 
        MG.Refine();
        Stokes.CreateNumberingVel(MG.GetLastLevel(), vidx1);    
        Stokes.CreateNumberingPr(MG.GetLastLevel(), pidx1);    
        std::cerr << "altes und neues TriangLevel: " << vidx2->TriangLevel << ", "
                  << vidx1->TriangLevel << std::endl;
        MG.SizeInfo(std::cerr);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        std::cerr << "Anzahl der Druck-Unbekannten: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cerr << "Anzahl der Geschwindigkeitsunbekannten: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;
        if (p2->RowIdx)
        {
            const StokesBndDataCL& BndData= Stokes.GetBndData();
            P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &BndData.Vel, &MG);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &BndData.Vel, &MG);
            Interpolate(pr1, pr2);
            Interpolate(vel1, vel2);
//            CheckSolution(v1,p1,&LsgVel,&LsgPr);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        time.Reset();
        time.Start();
        Stokes.SetupSystem(A, b, B, c);
        time.Stop();
        std::cerr << time.GetTime() << " seconds for setting up all systems!" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cerr << " A*x took " << time.GetTime() << " seconds!" << std::endl;
        time.Reset();
        time.Start();
        transp_mul( A->Data, v1->Data);
        time.Stop();
        std::cerr << "AT*x took " << time.GetTime() << " seconds!" << std::endl;
        
        { // write system in files for MatLab
            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
            Adat << A->Data;   Bdat << B->Data;    bdat << b->Data;    cdat << c->Data;
        }
        Stokes.GetDiscError(&LsgVel, &LsgPr);
//std::cout << A->Data << std::endl << b->Data << std::endl;
/*        double half= M_PI/8;
        MultiGridCL::TriangVertexIteratorCL vert= MG.GetTriangVertexBegin(A->RowIdx->TriangLevel);
        while (vert->GetCoord()[0]!=half || vert->GetCoord()[1]!=half || vert->GetCoord()[2]!=half) ++vert;
        IdxT unk= vert->Unknowns(A->RowIdx->Idx)[0];
        std::cerr << vert->GetCoord() << " has index " << unk << std::endl;
        std::cerr << "A(i,i) = " << A->Data(unk,unk) <<std::endl;    
        std::cerr << "B(i,j) = " << B->Data(vert->Unknowns(B->RowIdx->Idx)[0],unk) << std::endl;
*/        Uint meth;
//        int max_iter;
//        double tol;
        std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
        time.Reset();

        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        Stokes.SetupMass( &M);

        double outer_tol;
        std::cerr << "tol = "; std::cin >> outer_tol;

        if (meth)
        {
            PSchur_PCG_CL schurSolver( M.Data, outer_tol, 200, inner_iter_tol, 200);
            time.Start();
            schurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
        }
        else // Uzawa
        {
            double tau;
            Uint inner_iter;
            std::cerr << "tau = "; std::cin >> tau;
            std::cerr << "#PCG steps = "; std::cin >> inner_iter;
            Uzawa_PCG_CL uzawaSolver( M.Data, outer_tol, 5000, inner_iter_tol, inner_iter, tau);
            time.Start();
            uzawaSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "Iterationen: " << uzawaSolver.GetIter()
                      << "Norm des Res.: " << uzawaSolver.GetResid() << std::endl;
        }
/*        if (meth)
        {
            SsorPcCL<VectorCL, double> pc(omega);
            MatDescCL M;
            M.SetIdx( pidx1, pidx1);
            Stokes.SetupMass( &M);
            SsorMassPcCL<SchurComplMatrixCL> schur_pc( M.Data);
            SchurComplMatrixCL BABT(A->Data, B->Data, inner_iter_tol, omega);
            double outer_tol;
            std::cerr << "tol = "; std::cin >> outer_tol;
            time.Start();
            VectorCL rhs= -c->Data;
            {
                double tol= inner_iter_tol;
                int max_iter= 200;        
                VectorCL tmp(vidx1->NumUnknowns);
                PCG(A->Data, tmp, b->Data, pc, max_iter, tol);
                std::cerr << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;
                rhs+= B->Data*tmp;
            }
            std::cerr << "rhs has been set!" << std::endl;
            max_iter= 200;   
            tol= outer_tol;     

            PCG(BABT, p1->Data, rhs, schur_pc, max_iter, tol);
            std::cerr << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;

            tol= outer_tol;
            max_iter= 200;        
            PCG(A->Data, v1->Data, b->Data - transp_mul(B->Data, p1->Data), pc, max_iter, tol);
            time.Stop();
        }
        else
        {
            max_iter= 5000;
            double tau;
            Uint inner_iter;
            std::cerr << "tol = "; std::cin >> tol;
            std::cerr << "tau = "; std::cin >> tau;
            std::cerr << "#PCG steps = "; std::cin >> inner_iter;
            MatDescCL I;
            I.SetIdx( pidx1, pidx1);
            Stokes.SetupMass( &I);
            time.Start();
            Uzawa( A->Data, B->Data, I.Data, v1->Data, p1->Data, b->Data, c->Data, tau, max_iter, tol, inner_iter, inner_iter_tol);
//            Uzawa2( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data, max_iter, tol, inner_iter, inner_iter_tol);
            time.Stop();
            std::cerr << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;
        }
*/        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
        Stokes.CheckSolution(v1, p1, &LsgVel, &LsgPr);
        A->Reset();
        B->Reset();
        b->Reset();
        c->Reset();
//        std::cerr << "Loesung Druck: " << p1->Data << std::endl;
//        CreateNumbering(new_idx->TriangLevel, err_idx);
//        err->SetIdx(err_idx);
//        NumMarked= EstimateError(new_x, 0.2, &Estimator);
//TODO: Fehler schaetzen
//        new_marks= EstimateError(new_x, 1.5, akt_glob_err, &ResidualErrEstimator);
//        err->Reset();
//        DeleteNumbering(err_idx);
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cerr << std::endl;
    }
    while (++step<maxStep);
    // we want the solution to be in _v
    if (v2 == &loc_v)
    {
        Stokes.v.SetIdx(&Stokes.vel_idx);
        Stokes.p.SetIdx(&Stokes.pr_idx);
        *Stokes.v.RowIdx= *loc_v.RowIdx;
        *Stokes.p.RowIdx= *loc_p.RowIdx;
        
        Stokes.v.Data.resize(loc_v.Data.size());
        Stokes.v.Data= loc_v.Data;
        Stokes.p.Data.resize(loc_p.Data.size());
        Stokes.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

int main (int argc, char** argv)
{
  try
  {
    if (argc!=3)
    {
        std::cerr << "You have to specify two parameters:\n\tsdrops <omega> <inner_iter_tol>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= M_PI/4.;

    typedef DROPS::StokesP2P1CL<DROPS::BrickBuilderCL, StokesCoeffCL> 
            StokesOnBrickCL;
    typedef StokesOnBrickCL MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel, &LsgVel};
        
    StokesOnBrickCL prob(brick, StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;
//    MarkAll(mg);


//    mg.Refine();
//   MarkAll(mg);
//    mg.Refine();
/*
    for (DROPS::Uint i=0; i<4; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        MarkAll(mg);
//        MarkSome(mg);
        std::cerr << DROPS::SanityMGOutCL(mg);
//        std::cerr << DROPS::DumpMGCL(mg);
//        std::cerr << i << std::endl;
        mg.Refine();
    }
*/
//    MarkAll(mg);
    double omega= atof(argv[1]);
    double inner_iter_tol= atof(argv[2]);
    std::cerr << "Omega: " << omega << " inner iter tol: " << inner_iter_tol << std::endl;
    Strategy(prob, omega, inner_iter_tol, 4);
    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    fil << DROPS::GeomSolOutCL<MyStokesCL::DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;

//    std::cout.flags(std::ios::fixed|std::ios::showpoint);
//    for ( DROPS::MultiGridCL::const_TriangVertexIteratorCL tit=const_cast<const DROPS::MultiGridCL&>(mg).GetTriangVertexBegin(prob.GetSolution().GetTriangLevel()),
//          theend=const_cast<const DROPS::MultiGridCL&>(mg).GetTriangVertexEnd(prob.GetSolution().GetTriangLevel()); tit!=theend; ++tit )
//        std::cout << prob.GetSolution().GetVertSol(*tit) << '\n';
//    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
//    std::cerr << "Verts: "   << mg.GetVertices().GetFullSize()
//              << " Edges: "  << mg.GetEdges().GetFullSize()
//              << " Tetras: " << mg.GetTetras().GetFullSize()
//              << std::endl;
/*
//int wait;
//cin>>wait;
    for (DROPS::Uint i=0; i<6; ++i)
    {
        UnMarkDrop(mg, min(mg.GetLastLevel(),10u));
//        std::cerr << DROPS::SanityMGOutCL(mg);
//        std::cerr << i << std::endl;
        mg.Refine();
    }
    UnMarkAll(mg);
    mg.Refine();
    UnMarkAll(mg);
    mg.Refine();
    UnMarkAll(mg);
    mg.Refine();
    UnMarkAll(mg);
    mg.Refine();
*/
//    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
/*
    std::cerr << "Verts: "   << mg.GetVertices().GetFullSize()
         << " Edges: "  << mg.GetEdges().GetFullSize()
         << " Tetras: " << mg.GetTetras().GetFullSize()
         << std::endl;
*/
//    std::cout << DROPS::GeomMGOutCL(mg, -1, true) << std::endl;
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

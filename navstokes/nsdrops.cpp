#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "navstokes/navstokes.h"
#include <fstream>


struct NS1CL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]= p[0];
        ret[1]= p[1];
        ret[2]= -2*p[2];
        return ret;
    }

    static double LsgPr(const DROPS::Point3DCL& p)
    {
        return (p[0]*p[0]+p[1]*p[1]+p[2]*p[2])/2;
    }

    // boundary value functions (in 2D-bnd-coords)
    //DROPS::StokesVelBndDataCL::bnd_type bnd_val_e2e3(const Point2DCL&);
    //DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e3(const Point2DCL&);
    //DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e2(const Point2DCL&);


    // q*u - nu*laplace u + (u*D)u - Dp = f
    //                            div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p)
            { DROPS::SVectorCL<3> ret(0.0); ret[2]= 3*p[2]; return ret; }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };

};

struct NS2CL
{
    static DROPS::SVectorCL<3> LsgVel(const DROPS::Point3DCL& p)
    {
        DROPS::SVectorCL<3> ret;
        ret[0]= 1;
        ret[1]= 1;
        ret[2]= 1;
        return ret;
    }

    static double LsgPr(const DROPS::Point3DCL& p)
    {
        return 0;
    }

    // boundary value functions (in 2D-bnd-coords)
    //DROPS::StokesVelBndDataCL::bnd_type bnd_val_e2e3(const Point2DCL&);
    //DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e3(const Point2DCL&);
    //DROPS::StokesVelBndDataCL::bnd_type bnd_val_e1e2(const Point2DCL&);


    // q*u - nu*laplace u + (u*D)u - Dp = f
    //                            div u = 0
    class StokesCoeffCL
    {
      public:
        static double q(const DROPS::Point3DCL&) { return 0.0; }
        static DROPS::SVectorCL<3> f(const DROPS::Point3DCL& p)
            { DROPS::SVectorCL<3> ret(0.0); return ret; }
        const double nu;

        StokesCoeffCL() : nu(1.0) {}
    };

};

typedef NS1CL MyPdeCL;

namespace DROPS // for Strategy
{

template<class PreCondT, class SchurPreCondT>
void Schur( const MatrixCL& A, const PreCondT& pc, const MatrixCL& B, 
            VectorCL& u, VectorCL& p, const VectorCL& b, const VectorCL& c,
            const SchurPreCondT& schur_pc,
            const double inner_tol, const double outer_tol, const Uint max_iter, const Uint outer_iter)
// solve:       S*q = B*(A^-1)*b - c
//              A*u = b - BT*q
//                q = dt*p 
{
    VectorCL rhs= -c;
    {
        double tol= inner_tol;
        int iter= max_iter;        
        VectorCL tmp( u.size());
        PCG( A, tmp, b, pc, iter, tol);
        std::cerr << "Iterationen: " << iter << "    Norm des Residuums: " << tol << std::endl;
        rhs+= B*tmp;
    }
    std::cerr << "rhs has been set! Now solving pressure..." << std::endl;
    int iter= outer_iter;   
    double tol= outer_tol;     
    PCG( SchurComplMatrixCL( A, B, inner_tol, 1), p, rhs, schur_pc, iter, tol);
    std::cerr << "Iterationen: " << iter << "    Norm des Residuums: " << tol << std::endl;

    std::cerr << "pressure has been solved! Now solving velocities..." << std::endl;
    tol= outer_tol;
    iter= max_iter;        
    PCG(A, u, b - transp_mul(B, p), pc, iter, tol);
    std::cerr << "Iterationen: " << iter << "    Norm des Residuums: " << tol << std::endl;
    std::cerr << "-----------------------------------------------------" << std::endl;
}

template<class MGB, class Coeff>
void Strategy(NavierStokesP2P1CL<MGB,Coeff>& NS, double inner_iter_tol, Uint maxStep)
// flow control
{
    typedef typename NavierStokesP2P1CL<MGB,Coeff>::VelVecDescCL VelVecDescCL;
    
    MultiGridCL& MG= NS.GetMG();

    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= &NS.vel_idx;
    IdxDescCL* pidx1= &NS.pr_idx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx2= &loc_pidx;
//    IdxDescCL* err_idx= &_err_idx;
    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &NS.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &NS.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &NS.b;
    VelVecDescCL* c= &NS.c;
//    VecDescCL* err= &_err;
    MatDescCL* A= &NS.A;
    MatDescCL* B= &NS.B;
    MatDescCL* N= &NS.N;
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
        NS.CreateNumberingVel(MG.GetLastLevel(), vidx1);    
        NS.CreateNumberingPr(MG.GetLastLevel(), pidx1);    
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
        if (v2->RowIdx)
        {
            const StokesBndDataCL& BndData= NS.GetBndData();
            P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
            P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &BndData.Vel, &MG);
            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &BndData.Vel, &MG);
            Interpolate(pr1, pr2);
            Interpolate(vel1, vel2);
//            CheckSolution(v1,p1,&MyPdeCL::LsgVel,&MyPdeCL::LsgPr);
            v2->Reset();
            p2->Reset();
        }
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        N->SetIdx(vidx1, vidx1);
        time.Reset();
        time.Start();
        NS.SetupSystem(A, b, B, c);
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
/*        
        { // write system in files for MatLab
            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), Ndat("Nmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
            Adat << A->Data;   Bdat << B->Data; Ndat << N->Data;   bdat << b->Data;    cdat << c->Data;
        }
*/      NS.GetDiscError(&MyPdeCL::LsgVel, &MyPdeCL::LsgPr);
        Uint meth, schritt= 0;
        int numsteps, outer_iter;
        double tol, outer_tol;
        std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
        std::cerr << "tol = "; std::cin >> tol;
        if (!meth) // Uzawa
            { std::cerr << "# PCG steps = "; std::cin >> numsteps; }
        else
            { std::cerr << "# outer CG steps = "; std::cin >> numsteps; }
        
        // adaptive fixedpoint defect correction
        //---------------------------------------
        time.Reset();
        VectorCL d( vidx1->NumUnknowns), e( pidx1->NumUnknowns),
                 w( vidx1->NumUnknowns), q( pidx1->NumUnknowns);
        VelVecDescCL rhsN( vidx1), v_omw( vidx1);
        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        NS.SetupMass( &M);
        double omega= 1, res; // initial value (no damping)
        SsorPcCL<VectorCL, double> pc(1.);
        SsorMassPcCL<SchurComplMatrixCL> schur_pc( M.Data);
        for(;;) // ever
        {
            NS.SetupNonlinear( N, v1, &rhsN);
            MatrixCL AN;
            ConvexComb( AN, 1., A->Data, 1., N->Data);
            
            // calculate defect:
            d= AN*v1->Data + transp_mul( B->Data, p1->Data) - b->Data - rhsN.Data;
            e= B->Data*v1->Data                             - c->Data;
            
            std::cerr << (++schritt) << ": res = " << (res= sqrt(d*d + e*e) ) << std::endl; 
            if (res < tol )
                break;
            
            // solve correction:
            outer_iter= meth ? numsteps : 500;
            outer_tol= res/1000;
            if (outer_tol < tol) outer_tol= tol;
            if (meth)
                Schur( AN, pc, B->Data, w, q, d, e, schur_pc, inner_iter_tol, outer_tol, 100, outer_iter);
            else
            {
                Uzawa( AN, B->Data, M.Data, w, q, d, e, 1, outer_iter, outer_tol, numsteps, inner_iter_tol);
                std::cerr << "iteration stopped after step " << outer_iter 
                          << " with res = " << outer_tol << std::endl;
            }
            
            // calculate adaption:
            N->Data.clear();
            v_omw.Data= v1->Data - omega*w;
            NS.SetupNonlinear( N, &v_omw, &rhsN);
            
            d= A->Data*w + N->Data*w + transp_mul( B->Data, q);
            e= B->Data*w;
            omega= d*(A->Data*v1->Data) + d*(N->Data*v1->Data) + d*transp_mul( B->Data, p1->Data) + e*(B->Data*v1->Data)
                 - d*b->Data - d*rhsN.Data - e*c->Data;
            omega/= d*d + e*e;
            std::cerr << "omega = " << omega << std::endl;
            
            // update solution:
            v1->Data-= omega*w;
            p1->Data-= omega*q;
        }
        time.Stop();
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
        
/*        time.Reset();
        if (meth)
        {
            SsorPcCL<VectorCL, double> pc(omega);
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
//            tol= 1.0e-14;
            max_iter= 200;   
            tol= outer_tol;     
    //        PCG(A->Data, new_x->Data, b->Data, pc, max_iter, tol);
            CG(BABT, p1->Data, rhs, max_iter, tol);
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
            NS.SetupMass( &I);
            time.Start();
            Uzawa( A->Data, B->Data, I.Data, v1->Data, p1->Data, b->Data, c->Data, tau, max_iter, tol, inner_iter, inner_iter_tol);
//            Uzawa2( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data, max_iter, tol, inner_iter, inner_iter_tol);
            time.Stop();
            std::cerr << "Iterationen: " << max_iter << "    Norm des Residuums: " << tol << std::endl;
        }
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
*/        NS.CheckSolution(v1, p1, &MyPdeCL::LsgVel, &MyPdeCL::LsgPr);
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
        NS.v.SetIdx(&NS.vel_idx);
        NS.p.SetIdx(&NS.pr_idx);
        *NS.v.RowIdx= *loc_v.RowIdx;
        *NS.p.RowIdx= *loc_p.RowIdx;
        
        NS.v.Data.resize(loc_v.Data.size());
        NS.v.Data= loc_v.Data;
        NS.p.Data.resize(loc_p.Data.size());
        NS.p.Data= loc_p.Data;
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
    if (argc!=2)
    {
        std::cerr << "You have to specify one parameter:\n\tnsdrops <inner_iter_tol>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= M_PI/4.;

    typedef DROPS::NavierStokesP2P1CL<DROPS::BrickBuilderCL, MyPdeCL::StokesCoeffCL> 
            NSOnBrickCL;
    typedef NSOnBrickCL MyNavierStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 1, 1, 1);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel, &MyPdeCL::LsgVel};
        
    MyNavierStokesCL prob(brick, MyPdeCL::StokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;
    double inner_iter_tol= atof(argv[1]);
    std::cerr << "inner iter tol: " << inner_iter_tol << std::endl;
    
    Strategy(prob, inner_iter_tol, 4);

    std::cerr << "hallo" << std::endl;
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    fil << DROPS::GeomSolOutCL<MyNavierStokesCL::DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, min, max) << std::endl;

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

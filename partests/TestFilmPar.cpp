//**************************************************************************
// File:    TestFilmPar.cpp                                                *
// Content: parallel solver for a film problem                             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   04. July 2007                                                  *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestFilmPar.cpp
/// \brief parallel solver for a film problem

 // include parallel computing!
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/exchange.h"          // transfer of numerical datas
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construuct the initial multigrid
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid

 // include numeric computing!
#include "num/parsolver.h"              // various parallel solvers
#include "num/parprecond.h"             // various parallel preconditioners
#include "num/parlanczos.h"
#include "num/fe.h"

 // include problem class
#include "poisson/instatpoisson.h"      // setting up the Poisson problem

 // include misc
#include "out/output.h"                 // output in geomview or other formats
#include "out/ensightOut.h"             // output in ensight format
#include "partests/params.h"            // reading parameter-files

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;
/****************************************************************************
* G L O B A L  V A R I B L E S                                              *
****************************************************************************/
DROPS::TimeStoreCL Times(2);   // time table all times are stored within this class
enum TimePart{                  // numbers for accesing the time table
    T_Discretize,
    T_Solve
};

/****************************************************************************
* S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_Discretize, "Discretize");
    Times.SetDescriber(T_Solve, "Solve");
}

DROPS::ParamParFilmCL C;     // parameter class
const char line[] ="----------------------------------------------------------------------------------";
using DROPS::ProcCL;

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
****************************************************************************/
bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    ofstream check(dat);

    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);

    check.close();
    bool sane=DROPS::Check(pmg_sane && mg_sane);
    if (!sane)
        throw DROPS::DROPSErrCL("CheckParMultiGrid: Error within multigrid!");
    return sane;
}

/****************************************************************************
* P O I S S O N  C O E F F  C L                                             *
*****************************************************************************
*   Coeffs of the Poisson problem:                                          *
*       - Laplace u + Vel.(grad u) +q*u = f                                 *
*       with q=0 and Vel=0                                                  *
****************************************************************************/
class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&, double){
        return 0.0;
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL p, double)
    {
        const double g=C.g[1],
                     dy=C.Dim[1],
                     Re=std::pow(dy,3)*g / (3.*std::pow(C.KineticViscosity,2)),
//                      Re=50.,
                     factor=12.*Re*C.KineticViscosity/std::pow(dy,3),
                     y=p[1];

        DROPS::Point3DCL v(0.0);
        v[0]= factor*(y*dy/2. - std::pow(y,2)/2.);
        return v;
    }
    static double f(const DROPS::Point3DCL&, double= 0.0){
        return 0;
    }
    static double alpha(const DROPS::Point3DCL&, double=0.0){
        return 1.;
    }
};

double Inflow(const DROPS::Point3DCL& p, double)
{
    const double y=p[1];
    if (y>C.Dim[1]*0.05)
        return C.InflowTemperature;
    else
        return C.WallTemperature-(C.WallTemperature-C.InflowTemperature)/(C.Dim[1]*0.05)*y;

}

double Wall(const DROPS::Point3DCL& , double)
{
    return C.WallTemperature;
}

double Null(const DROPS::Point3DCL&, double)
{
    return 0.;
}

namespace DROPS
{
/****************************************************************************
* S O L V E                                                                 *
*****************************************************************************
*   Embrace all functions solve a linear system. This function constructs   *
*   all solvers, that can be used, but uses only the wanted one             *
****************************************************************************/
template <typename Mat, typename Vec>
void Solve(const Mat &A, Vec &x, const Vec &b, ExchangeCL &ExX)
{
    // Measure the used time for solving the linear equation system
    ParTimerCL time;

    // Definition and initialisation of GMRES-Solver with preconditioning
    typedef ParJac0CL<ExchangeCL> PrecT;
    PrecT PCJac(ExX, C.Relax);
    ParPreGMResSolverCL<PrecT,ExchangeCL>
        P_J_GMRESSolver(C.Restart, C.Iter, C.Tol, ExX, PCJac, false, true, C.UseMGS);

    if (ProcCL::IamMaster())
        std::cout <<" Solve the linear system with GMRES("<<C.Restart<<"),"
                  <<" tol "<<C.Tol<<", max_iter "<<C.Iter<<std::endl;

    // Start time measurement
    time.Reset();
    time.Reset();

    // Solve the system
    P_J_GMRESSolver.Solve(A, x, b);

    time.Stop();
    if (ProcCL::IamMaster())
        std::cout <<" Solved the system with\n"
                  <<" - resid "<<P_J_GMRESSolver.GetResid()<<'\n'
                  <<" - steps "<<P_J_GMRESSolver.GetIter()<<'\n'<<line<<std::endl;
    Times.AddTime(T_Solve, time.GetMaxTime());
}


template<class PoissonCoeffCL>
void WriteEnsightFiles(const MultiGridCL& mg, const InstatPoissonP1CL<PoissonCoeffCL>& Poisson)
{
    EnsightP2SolOutCL ensight( mg, Poisson.idx.GetFinestPtr(), /*binary=*/false, /*masterout=*/true);

    const string filename= C.EnsDir + "/" + C.EnsCase;
    const string datgeo= filename+".geo",
                 dattmp= string(filename+".tmp");
    ensight.CaseBegin( string(C.EnsCase+".case").c_str());
    ensight.DescribeGeom(   C.geomName.c_str(), datgeo);
    ensight.DescribeScalar( C.varName.c_str(),  dattmp);
    ensight.Commit();
    ensight.CaseEnd();
    ensight.putGeom( datgeo);
    ensight.putScalar(dattmp, Poisson.GetSolution());
}

/****************************************************************************
* S T R A T E G Y                                                           *
*****************************************************************************
* The heattransfer problem is solved by this non adaptive procedure. The    *
* main steps are:                                                           *
*   1. Numbering the unknwons (map unknowns to geometry) and create         *
*      ExchangeCL (for communication of the unknowns)                       *
*   2. Discretize (setting up the linear equation system)                   *
*   3. Solve the linear equation system                                     *
****************************************************************************/
template<class PoissonCoeffCL>
void Strategy(InstatPoissonP1CL<PoissonCoeffCL>& Poisson)
{
    // Time meassurement
    ParTimerCL time;
    if (ProcCL::IamMaster()){
        std::cout <<line<< "\n Initialize problem, vectors and matrices " << std::endl;
    }

    // Problem class
    typedef InstatPoissonP1CL<PoissonCoeffCL> MyPoissonCL;

    // Geometry
    MultiGridCL& MG= Poisson.GetMG();

    // Mapping unknowns onto vertices
    MLIdxDescCL* idx= &Poisson.idx;

    // matrices and vectors of the linear equation system:
    VecDescCL* T= &Poisson.x;   VecDescCL* b= &Poisson.b;
    MLMatDescCL* A= &Poisson.A; MLMatDescCL* M= &Poisson.M;

    // helper matrices and vectors
    VecDescCL vU, vA, vM, vf;
    MLMatDescCL* U= &Poisson.U;
    MLMatrixCL AU, AUM;

    // Communication class
    ExchangeCL &ExX = Poisson.GetEx();

    // Point 1
    //--------
    // Create numbering of the unknowns according to P1 finite elements. The
    // ExchangeCL within the PoissonCL is also initialised.
    idx->SetFE( P1_FE);
    Poisson.CreateNumbering(MG.GetLastLevel(), idx);

    // tell vectors and matrices about the unknowns
    b->SetIdx( idx); T->SetIdx( idx); vU.SetIdx( idx);
    vA.SetIdx( idx); vM.SetIdx( idx); vf.SetIdx( idx);
    A->SetIdx(idx, idx); M->SetIdx(idx, idx); U->SetIdx(idx, idx);

    IdxT numUnknowns=idx->GetGlobalNumUnknowns(MG);
    if (ProcCL::IamMaster())
        std::cout << " Number of unknowns: " << numUnknowns << std::endl;

    // Point 2
    //--------
    if (ProcCL::IamMaster()){
        std::cout << line << std::endl;
        std::cout << " Discretize ... " << std::endl;
    }

    time.Reset();
    time.Reset();
    Poisson.SetupInstatSystem(*A, *M,0.);
    Poisson.SetupConvection(*U,vU,0.);
    Poisson.SetupInstatRhs(vA,vM,0.,vf,0.);
    AU.LinComb( C.HeatConductivity, A->Data, 1., U->Data);
    AUM.LinComb( 1., AU, 1., M->Data);
    b->Data = vf.Data + vU.Data + vA.Data  + vM.Data;
    time.Stop(); ::Times.AddTime(T_Discretize, time.GetMaxTime());

    // Point 3
    //--------
    // Solve the linear equation system
    if (ProcCL::IamMaster())
        std::cout << line << std::endl;
    Solve(AUM, T->Data, b->Data, ExX);


    if (C.Ensight){
        if (ProcCL::IamMaster())
            std::cout << " Write ensight-Data files" << std::endl;
        WriteEnsightFiles(MG, Poisson);
    }
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
    DROPS::ProcCL Proc(&argc, &argv);
    try
    {
        SetDescriber();

        if (argc!=2 && ProcCL::IamMaster()){
            std::cout << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param && ProcCL::IamMaster()){
            std::cout << "error while opening parameter file\n"; return 1;
        }
        param >> C; param.close();
        if (ProcCL::IamMaster())
            std::cout << C << std::endl;

        DROPS::ParTimerCL time, alltime;

        typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL> PoissonOnBCL;
        typedef PoissonOnBCL                             MyPoissonCL;

        // Init of the parallel structurs. Tell the ParMultiGrid class how many indices should be handled.
        DROPS::ParMultiGridCL pmg(0);

        // origin and dimension of the brick
        DROPS::Point3DCL orig(0.);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]=C.Dim[0];
        e2[1]=C.Dim[1];
        e3[2]=C.Dim[2];


        // bc={left, right, down, up, front, back}
        const DROPS::BndCondT bc[6]=
            { DROPS::DirBC, DROPS::OutflowBC,
              DROPS::DirBC, DROPS::OutflowBC,
              DROPS::OutflowBC, DROPS::OutflowBC};

        const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
            { &Inflow, &Null,
              &Wall,   &Null,
              &Null,   &Null};

        DROPS::InstatPoissonBndDataCL bdata(6, bc, bnd_fun);

        if (ProcCL::IamMaster()){
            DROPS::Point3DCL y; y[1]=(C.Dim[1])/2.;
            std::cerr << "Maximale Geschwindigkeit: "<<(PoissonCoeffCL::Vel(y,0.))<<" eps="<<C.HeatConductivity<<std::endl;
        }

        if (ProcCL::IamMaster()){
            std::cout << line << '\n'
                 << " Create initial grid and distribution ... \n";
        }

        // Create the MultiGrid on a single processor
        DROPS::MGBuilderCL * mgb;
        if (ProcCL::IamMaster())
            mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, C.nx, C.ny, C.nz);
        else
            mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);

        // Setup the problem
        PoissonOnBCL prob(*mgb, PoissonCoeffCL(), bdata);
        DROPS::MultiGridCL &mg = prob.GetMG();
        pmg.AttachTo(mg);

        // Init the LoadBalHandler (create an MultiGridCL and distribute the multigrid)
        DROPS::LoadBalHandlerCL lb(mg);
        lb.DoInitDistribution(ProcCL::Master());
        for (int ref=0; ref<C.Refine; ++ref){
            DROPS::MarkAll(mg);
            pmg.Refine();
            lb.DoMigration();
        }

        if (ProcCL::IamMaster())
            cout << " Check the parallel multigrid ... " << flush;
        if (CheckParMultiGrid(pmg)){
            if (ProcCL::IamMaster()) cout << " OK\n";
        }
        else{
            if (ProcCL::IamMaster()) cout << "  !!! not OK !!!\n";
            throw DROPS::DROPSErrCL("Error in parallel multigrid");
        }

        if (ProcCL::IamMaster())
            cout <<line<< "\n Distribution of elemente:\n";
        mg.SizeInfo(cout);


        DROPS::Strategy(prob);

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        Times.Print(cout);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}

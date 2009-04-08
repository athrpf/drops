//**************************************************************************
// File:    TestGMRESSer.cpp                                               *
// Content: Sequential solver for Poisson problems                         *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   10. December 2007                                              *
//**************************************************************************
// test different setting and solvers for different multigrids for         *
// parallel solving of Poisson problems                                    *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestPoissonSer.cpp
/// \brief Testing sequential solving of the Poisson problem


 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construuct the initial multigrid

 // include numeric computing!
#include "num/solver.h"                 // various sequential solvers
#include "num/fe.h"

 // include problem class
#include "poisson/instatpoisson.h"      // setting up the Poisson problem

 // include misc
#include "out/output.h"                 // output in geomview or other formats
#include "partests/params.h"            // reading parameter-files

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

// parameter settings for the Poisson problem and solver class
DROPS::ParamParPoissonCL C;
const char line[] ="----------------------------------------------------------------------------------";

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
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL, double){
        DROPS::Point3DCL ret; //ret[0]=1.; ret[1]=1.; ret[2]=1.;
        return ret;
    }
    static double f(const DROPS::Point3DCL& p, double= 0.0){
        return 128*(p[1]*p[2]*(1.-p[1])*(1.-p[2])
                + p[0]*p[2]*(1.-p[0])*(1.-p[2])
                + p[0]*p[1]*(1.-p[0])*(1.-p[1]));
    }
    static double alpha(const DROPS::Point3DCL&, double=0.0){
        return 1.;
    }
};

/****************************************************************************
* S O L U T I O N                                                           *
*****************************************************************************
*   Solution of the poisson problem to the above coeffs with Dirichlet      *
*   boundary (zero)                                                         *
****************************************************************************/
inline double Lsg( const DROPS::Point3DCL& p, double=0.0){
    return 64*p[0]*p[1]*p[2]*(1.-p[0])*(1.-p[1])*(1.-p[2]);
}

// boundary function
double Null (const DROPS::Point3DCL&, double) {return 0.;}

namespace DROPS
{

/****************************************************************************
* S T R A T E G Y                                                           *
*****************************************************************************
* The Poisson problem is solved by this non adaptive procedure. These are   *
* the main steps:                                                           *
*   1. Numbering the unknwons (map unknowns to geometry)                    *
*   2. Discretize                                                           *
*   3. Solve the linear equation system                                     *
****************************************************************************/
template<class PoissonCoeffCL>
void Strategy(InstatPoissonP1CL<PoissonCoeffCL>& Poisson)
{
    TimerCL time;

    typedef InstatPoissonP1CL<PoissonCoeffCL> MyPoissonCL;

    MultiGridCL& MG= Poisson.GetMG();

    IdxDescCL* idx= &Poisson.idx;
    VecDescCL* x= &Poisson.x; VecDescCL* b= &Poisson.b;
    VecDescCL vU, vA, vM, vf;
    MatDescCL* A= &Poisson.A; MatDescCL* M= &Poisson.M;
    MatDescCL* U= &Poisson.U;
    MatrixCL AU, AUM;

    // Point 1
    //--------
    idx->Set(1,0,0,0); // 1 unknown on vertex

    std::cout << line << std::endl;
    std::cout << " Numbering unknowns ... " << std::endl;

    // Create numbering of the unknowns according to the index.
    time.Reset();
    Poisson.CreateNumbering(MG.GetLastLevel(), idx);
    time.Stop();
    std::cout << " - Time for numbering unknowns "<<time.GetTime() << std::endl;
    std::cout << "   Number of unknwons: "<<idx->NumUnknowns << std::endl;

    // tell vectors and matrices about the unknowns
    b->SetIdx( idx); x->SetIdx( idx); vU.SetIdx( idx);
    vA.SetIdx( idx); vM.SetIdx( idx); vf.SetIdx( idx);
    A->SetIdx(idx, idx); M->SetIdx(idx, idx); U->SetIdx(idx, idx);

    // Point 2
    //--------
    std::cout << line << std::endl;
    std::cout << " Discretize ... " << std::endl;

    time.Reset();
    Poisson.SetupInstatSystem(*A, *M,0.);
    Poisson.SetupConvection(*U,vU,0.);
    Poisson.SetupInstatRhs(vA,vM,0.,vf,0.);
    AU.LinComb( C.nu, A->Data, 1., U->Data);
    AUM.LinComb( 1., AU, 1., M->Data);
    b->Data = vf.Data + vU.Data + vA.Data  + vM.Data;
    time.Stop();

    std::cout << " - Time for setting up matrices and rhs "<<time.GetTime()
              << " sec " << std::endl;

    // Point 3
    //--------
    // Solve
    std::cout << line << std::endl
              << " Solve the system ... " << std::endl;

    // typedefinition of the preconditioner
    typedef JACPcCL PreCondT;                       // Jacobi as preconditioner
    // typedefinition of solver
    typedef GMResSolverCL<PreCondT> SolverT;        // GMRES with preconditioning

    PreCondT pc;
    SolverT gmres(pc, C.restart, C.iter, C.tol, false, false, LeftPreconditioning);

    time.Reset();
    gmres.Solve(AUM, x->Data, b->Data);
    time.Stop();
    std::cout << " - Time for solving the system with GMRES " << time.GetTime() << '\n'
              << "   iter:  " << gmres.GetIter() << '\n'
              << "   resid: " << gmres.GetResid() << std::endl;
}

} // end of namespace DROPS

int main (int argc, char** argv)
{
    try
    {
        if (argc!=2){
            std::cout << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param){
            std::cout << "error while opening parameter file\n"; return 1;
        }
        param >> C; param.close();
        std::cout << C << std::endl;

        DROPS::TimerCL time;

        typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL> PoissonOnBCL;
        typedef PoissonOnBCL                             MyPoissonCL;

        DROPS::Point3DCL orig(0.);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= e2[1]= e3[2]= 1.0;

        // just Dirichlet boundarys
        const bool IsNeumann[6]=
            {false, false, false, false, false, false};
        const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
            { &Null, &Null, &Null, &Null, &Null, &Null};
        DROPS::InstatPoissonBndDataCL bdata(6, IsNeumann, bnd_fun);

        std::cout << line << std::endl
                  << " Create initial grid and set up problem ... \n";

        time.Reset();
        // Create discretized grid
        DROPS::BrickBuilderCL mgb(orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);

        // Setup the problem
        PoissonOnBCL prob(mgb, PoissonCoeffCL(), bdata);
        DROPS::MultiGridCL &mg = prob.GetMG();
        time.Stop();
        std::cout << " - Time for setting up problem "<<time.GetTime()<<" sec."<<std::endl;

        // Refine grid regular
        std::cout << line << std::endl
                  << " Refining grid ... " << std::endl;
        time.Reset();
        for (int ref=0; ref<C.refall; ++ref)
        {
            // Markieren und verfeinern
            std::cout << " Refine (" << ref << ") regular ... \n";
            DROPS::MarkAll(mg);
            mg.Refine();
        }
        time.Stop();
        std::cout << " - Time for refining grid took "<<time.GetTime()<<" sec."<<std::endl
                  << "  Number of geometry objects: ";

        mg.SizeInfo(cout);

        DROPS::Strategy(prob);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}

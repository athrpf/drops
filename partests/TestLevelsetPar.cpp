//**************************************************************************
// File:    TestLevelsetPar.cpp                                            *
// Content: Test enviroment for parallel fastmarching                      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   25. Januar 2006                                                *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestLevelsetPar.cpp
/// \brief Test enviroment for parallel fastmarching

 // include parallel computing!
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include "parallel/exchange.h"          // transfer of numerical datas
#include "parallel/parfastmarch.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construuct the initial multigrid
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "parallel/exchange.h"
#include "out/ensightOut.h"

 // include parallel fastmarching
#include "levelset/fastmarch.h"
#include "num/discretize.h"
#include "misc/problem.h"
#include "num/bndData.h"
#include "num/fe.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;
/****************************************************************************
* G L O B A L  V A R I B L E S                                              *
****************************************************************************/
DROPS::TimeStoreCL Times(3);    // time table all times are stored within this class
enum TimePart{                  // numbers for accesing the time table
    T_Init,
    T_Reparam,
    T_ReparamEuklid
};
const char line[] ="----------------------------------------------------------------------------------\n";

/****************************************************************************
* S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_Init, "Initialization");
    Times.SetDescriber(T_Reparam, "Reparametrization");
    Times.SetDescriber(T_ReparamEuklid, "Eukilidian Reparametrization");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/****************************************************************************
* D E F I N I T I O N   O F   T H E   P R O B L E M                         *
*****************************************************************************
* Only natural homogenious boundray conditions, no coeffs                   *
****************************************************************************/
class EmptyCoeffCL {} EmptyCoeff;
typedef DROPS::BndDataCL<double>                                        BndDataT;
typedef DROPS::ProblemCL<EmptyCoeffCL, BndDataT>                        TestP2LevelSetCL;
typedef DROPS::P2EvalCL<double, const BndDataT, const DROPS::VecDescCL> const_DiscSol;

using DROPS::ProcCL;

const int MIG=0, REF=1;         // constants for choosing the the filename in PrintMG
/****************************************************************************
* P R I N T  M U L T I G R I D  I N  A S C I I  F I L E                     *
****************************************************************************/
void PrintMG(const DROPS::ParMultiGridCL& pmg, int type=MIG)
{
    const int me=DROPS::ProcCL::MyRank();
    static int REFnum=0;
    static int MIGnum=0;
    char filename[30];

    if (type==REF)
        sprintf(filename, "output/%i_MG_REF_%i.mg",me,REFnum++);
    else if (type == MIG)
        sprintf(filename, "output/%i_MG_MIG_%i.mg",me,MIGnum++);
    if (me==0)
        std::cout << " - Writing multigrid into: " << filename<< " (for proc 0)"<<std::endl;

    ofstream file(filename);
    pmg.DebugInfo(file);
    file.close();
}

/****************************************************************************
* P R I N T   S O L U T I O N   F O R   E N S I G H T                       *
****************************************************************************/
// void PrintEnsight(const TestP2LevelSetCL& prob, const DROPS::VecDescCL& phi,
//                   const DROPS::VecDescCL& phiEukild,
//                   const DROPS::VecDescCL& phiDiff,
//                  )
// {
//     DROPS::EnsightP2SolOutCL ensight( prob.GetMG(), phi.RowIdx, false);
//
//     const std::string EnsCase ="Levelset",
//                       filename= "ensight/" + EnsCase,
//                       datgeo  = filename+".geo",
//                       dattmp  = (filename+".lst"),
//                       geomName= "Wuerfel",
//                       varName1 = "FastMarch",
//                       varName2 = "euklid",
//                       varName3 = "difference",
//
//     ensight.CaseBegin( string(EnsCase+".case").c_str());
//     ensight.DescribeGeom(   geomName.c_str(), datgeo);
//     ensight.DescribeScalar( varName.c_str(), dattmp);
//     ensight.Commit();
//     ensight.CaseEnd();
//     ensight.putGeom( datgeo);
//     const_DiscSol sol (&phi, &prob.GetBndData(), &prob.GetMG());
//     ensight.putScalar(dattmp, sol);
// }

/****************************************************************************
* M A R K  D R O P                                                          *
*****************************************************************************
*  Mark tetras around the phase-boundary                                    *
****************************************************************************/
void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

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
    return DROPS::Check(pmg_sane && mg_sane);
}

/****************************************************************************
* D I S T F U N C T I O N S                                                 *
*****************************************************************************
* Distance functions:                                                       *
*   1 Brick is divided into two parts. The upper and the lower part. The    *
*     border has a distance of 0.25 from the lower domain boundary          *
*   2 A drop in the middle of the brick (Radius 0.25)                       *
****************************************************************************/
const double Radius = 0.25;
double DistFunction_1(const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL middle;
    middle[0]=middle[1]=middle[2]=0.5;
    return (p-middle).norm()-Radius;
}

const double BorderDist = 0.25;
double DistFunction_2(const DROPS::Point3DCL& p)
{
    return p[2]-BorderDist;
}

/****************************************************************************
* I N I T   L E V E L S E T                                                 *
*****************************************************************************
* Init the values of the function Phi according to a scalar valued function *
****************************************************************************/
void Init( DROPS::scalar_fun_ptr phi0, const DROPS::MultiGridCL& mg, DROPS::VecDescCL& Phi)
{
    const DROPS::Uint lvl= Phi.GetLevel(),
                      idx= Phi.RowIdx->GetIdx();

    for (DROPS::MultiGridCL::const_TriangVertexIteratorCL it= mg.GetTriangVertexBegin(lvl),
         end= mg.GetTriangVertexEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( it->GetCoord());
    }
    for (DROPS::MultiGridCL::const_TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(lvl),
         end= mg.GetTriangEdgeEnd(lvl); it!=end; ++it)
    {
        Phi.Data[it->Unknowns(idx)]= phi0( DROPS::GetBaryCenter( *it));
    }
}

void DisturbLevelSet(DROPS::VectorCL& phi)
{
    for (DROPS::Uint i=0; i<phi.size(); ++i){
        phi[i] = phi[i]*100.;
    }
}

namespace DROPS
{
void Strategy(TestP2LevelSetCL &prob)
{
    // Index Describer für je eine Unbekannte auf Vertices und Edges
    MLIdxDescCL idx( P1_FE);
    ParTimerCL time;

    bool euklid=true, fastmarch=true;

    // Erstellen einer Nummerierung und fuellen der ExchangeCL
    MultiGridCL& mg= prob.GetMG();
    idx.CreateNumbering(mg.GetLastLevel(), mg, prob.GetBndData());
    ExchangeCL&  ex= idx.GetEx();

    if (ProcCL::IamMaster())
        std::cout << " Create ExchangeCL mit Abbildung der Indices\n";

    ex.SizeInfo(std::cout);

    if (ProcCL::IamMaster()){
        std::cout << " Meine Nachbar-Prozesse lauten: ";
        ExchangeCL::ProcNumCT neigh=ex.GetNeighbors();
        for (ExchangeCL::ProcNumCT::const_iterator it(neigh.begin()), end(neigh.end()); it!=end; ++it)
            std::cout << *it << " ";
        std::cout << std::endl;
    }

    if (ProcCL::IamMaster())
        std::cout << line << " Erstellen der Levelsetfunktion\n";
    // Erstellen der Levelsetfunktion
    VecDescCL  phi(&idx);
    Init(DistFunction_1, mg, phi);


    bool ensightOut=true;
    DROPS::EnsightP2SolOutCL *ensight=0;
    const std::string EnsCase ="Levelset",
                    filename= "ensight/" + EnsCase,
                    datgeo  = filename+".geo",
                    dattmp1  = (filename+".fm"),
                    dattmp2  = (filename+".eukl"),
                    dattmp3  = (filename+".diff"),
                    geomName= "Wuerfel",
                    varName1 = "FastMarch",
                    varName2 = "euklid",
                    varName3 = "difference";

    if (ensightOut)
    {

        if (ProcCL::IamMaster())
            std::cout << "Create ensight stuff\n";

        ensight = new DROPS::EnsightP2SolOutCL( prob.GetMG(), phi.RowIdx, false);

        ensight->CaseBegin( std::string(EnsCase+".case").c_str());
        ensight->DescribeGeom(   geomName.c_str(), datgeo);
        ensight->DescribeScalar( varName1.c_str(), dattmp1);
        ensight->DescribeScalar( varName2.c_str(), dattmp2);
        ensight->DescribeScalar( varName3.c_str(), dattmp3);
        ensight->Commit();
        ensight->CaseEnd();
        ensight->putGeom( datgeo);
    }


    if (ProcCL::IamMaster())
        std::cout << " Stoere die Levelsetfunktion\n";
    DisturbLevelSet(phi.Data);

    VectorCL phiFastMarch(phi.Data.size());

    if (fastmarch){
        if (ProcCL::IamMaster())
            std::cout << " Reparametrisieren der Levelsetfunktion mittels serial fastmarching\n";
        FastMarchCL FastMarch(mg, phi, ex);

        time.Reset();
        for (int i=0; i<2; ++i){
            FastMarch.Reparam();
        }
        time.Stop();
        Times.AddTime(T_Reparam, time.GetMaxTime());

        if (ensightOut){
            const_DiscSol sol1 (&phi, &prob.GetBndData(), &prob.GetMG());
            ensight->putScalar(dattmp1, sol1);
        }


        phiFastMarch=phi.Data;

        if (ProcCL::IamMaster())
        std::cout << "Check, ob levelset Funktion akkumuliert ist\n";

        if ( ex.IsAcc(phi.Data) ){
            if (ProcCL::IamMaster())
                std::cout << "  -> Phi ist akkumuliert\n";
        }
        else{
            if (ProcCL::IamMaster())
                std::cout << "  -> Phi ist nicht akkumuliert\n";
        }
    }

    if (euklid){
        if (ProcCL::IamMaster())
            std::cout << " Stoere die Levelsetfunktion\n";
        DisturbLevelSet(phi.Data);

        FastMarchCL FastMarch(mg, phi, ex);
        if (ProcCL::IamMaster())
            std::cout << " Reparametrisieren der Levelsetfunktion mittels eukilidischem Abstand\n";
        time.Reset();
        FastMarch.ReparamEuklid();
        time.Stop();
        Times.AddTime(T_ReparamEuklid, time.GetMaxTime());

        if (ensightOut){
            if (ProcCL::IamMaster())
                std::cout << " Ausgabe nach ensight\n";

            const_DiscSol sol2 (&phi, &prob.GetBndData(), &prob.GetMG());
            ensight->putScalar(dattmp2, sol2);
        }

        if (fastmarch){
            phi.Data -= phiFastMarch;
            // Ensight Ausgabe
            if (ensightOut){
                const_DiscSol sol3 (&phi, &prob.GetBndData(), &prob.GetMG());
                ensight->putScalar(dattmp3, sol3);
            }
        }

        if (ProcCL::IamMaster())
            std::cout << "Check, ob levelset Funktion akkumuliert ist\n";

        if ( ex.IsAcc(phi.Data) ){
            if (ProcCL::IamMaster())
                std::cout << "  -> Phi ist akkumuliert\n";
        }
        else{
            if (ProcCL::IamMaster())
                std::cout << "  -> Phi ist nicht akkumuliert\n";
        }
    }
    if (ensight) delete ensight;
}
} // end of namespace DROPS

int main (int argc, char** argv)
{
    DROPS::ProcCL Proc(&argc, &argv);
    try
    {
        SetDescriber();
        DROPS::ParTimerCL alltime;

        if (argc!=2){
            if (ProcCL::IamMaster()){
                std::cerr << "usage: "<<argv[0]<<" <num_ref>" << std::endl;
            }
            throw DROPS::DROPSErrCL("No enough parameters are given!");
        }
        const int num_ref = std::atoi(argv[1]);

        // Init of the parallel structurs.
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

        DROPS::Point3DCL orig(0.);
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]= e2[1]= e3[2]= 1.0;

        if (ProcCL::IamMaster())
            std::cout << line << " Create initial grid and distribution ... \n";

        DROPS::MGBuilderCL * mgb;
        if (DROPS::ProcCL::IamMaster())
            mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, 4, 4, 4);
        else
            mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);

        // Setup boundary data
        DROPS::BndCondT bndCond[6] =
//             {DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC};
                {DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC,DROPS::Nat0BC};
        BndDataT Boundary(6, bndCond);

        // Setup the problem
        TestP2LevelSetCL prob(*mgb, EmptyCoeff, Boundary);

        // Get MultiGrid, attach to pmg, distribute ...
        DROPS::MultiGridCL &mg = prob.GetMG();
        pmg.AttachTo(mg);

        // Init the LoadBalHandler (create an MultiGridCL and distribute the multigrid)
        DROPS::LoadBalHandlerCL lb(mg);
        lb.DoInitDistribution(ProcCL::Master());
        lb.SetStrategy(DROPS::Adaptive);

        // Refinement
        if (DROPS::ProcCL::IamMaster())
            std::cout << line << " Refine the grid "<<num_ref<<" times regular\n";

        for (int ref=0; ref<num_ref; ++ref)
        {
            // Markieren und verfeinern
            if (DROPS::ProcCL::IamMaster()) std::cout << " Refine (" << ref << ") regular ...\n";
//             MarkDrop(mg, mg.GetLastLevel());
            DROPS::MarkAll(mg);
            pmg.Refine();
            lb.DoMigration();
            Times.IncCounter(lb.GetMovedMultiNodes());
        }

        if (DROPS::ProcCL::IamMaster()) cout << " Check the parallel multigrid ... " << flush;
        if (CheckParMultiGrid(pmg)){
            if (DROPS::ProcCL::IamMaster()) cout << " OK\n";
        }
            else
                if (ProcCL::IamMaster()) cout << "  !!! not OK !!!\n";

        if (ProcCL::IamMaster())
            cout << " Distribution of elemente:\n";
        mg.SizeInfo(cout);

        DROPS::Strategy(prob);

        if (DROPS::ProcCL::IamMaster()){
            cout << line;
        }

        alltime.Stop();
        Times.SetOverall(alltime.GetMaxTime());
        Times.Print(cout);

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }
}

//**************************************************************************
// File:    TestMGSerPar.cpp                                               *
// Content: Test enviroment for serialization of a multigrid               *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, Timo Henrich, SC RWTH Aachen                 *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   15. November 2007                                              *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestMGSerPar.cpp
/// \brief Test enviroment for serialization of a multigrid

 // include parallel computing!
#include "parallel/parallel.h"          // proc handling, reduce operations, ...
#include "parallel/partime.h"           // parallel time-messurement
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"             // multigrid on each processor
#include "geom/builder.h"               // construuct the initial multigrid
#include "parallel/parmultigrid.h"      // handle multigrid over different procs
#include "parallel/loadbal.h"           // distribute multigrid
#include "parallel/parmgserialization.h"

#include "misc/problem.h"
#include "num/fe.h"

#include "partests/params.h"
#include "out/output.h"
#include "out/vtkOut.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

// Parameter-Klasse
DROPS::ParamParSerCL C;

const char line[]="---------------------------------------------------------------\n";
const double tolerance=DROPS::DoubleEpsC;

/****************************************************************************
* M A R K I N G   P R O C E D U R E S                                       *
****************************************************************************/
void MarkDrop (const DROPS::MultiGridCL& mg)
{
    DROPS::Point3DCL Mitte;
    Mitte[0]=0.5*(C.dim[0]-C.orig[0]);
    Mitte[1]=0.5*(C.dim[1]-C.orig[1]);;
    Mitte[2]=0.5*(C.dim[2]-C.orig[2]);;

    for (DROPS::MultiGridCL::const_TriangTetraIteratorCL It(mg.GetTriangTetraBegin()), ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
}

void MarkCorner (DROPS::MultiGridCL& mg)
{
    DROPS::Point3DCL Corner(C.orig);
    /// \todo Aufrufparameter maxLevel entfernt
    ///  - for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
    ///  - ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    ///waren die ursprünglichen Zeilen. So okay ?
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin()), ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
            if ( (GetBaryCenter(*It)-Corner).norm()<=0.3)
                It->SetRegRefMark();
}

/****************************************************************************
* F U N C T I O N S   C O N C E R N I N G   D O F S                         *
****************************************************************************/

/// Init pressure
double Pressure(const DROPS::Point3DCL& x, double)
{
    return (x-C.orig+C.dim).norm();
}

/// Init levelset as distance to origin of the brick
double Levelset(const DROPS::Point3DCL& x, double)
{
    return (x-C.orig).norm()+1.;
}

/// Init velocity as parabolic profile with zero velocities at the bottom of the brick
DROPS::Point3DCL Velocity(const DROPS::Point3DCL& x, double =0.)
{
    const double y=  (C.orig[0]-x[0])*(C.orig[0]+C.dim[0]-x[0])
                   * (C.orig[2]-x[2])*(C.orig[2]+C.dim[2]-x[2])
                   * (x[1]-C.orig[1]);
    return DROPS::MakePoint3D(10., y, 10.);
}

typedef double (*scal_fun)( const DROPS::Point3DCL&, double);
typedef DROPS::Point3DCL (*vec_fun)( const DROPS::Point3DCL&, double);

/// Set scalar valued function on vertices and edges
void SetDOF(const DROPS::MultiGridCL& mg, DROPS::VecDescCL& vec, scal_fun f)
{
    using namespace DROPS;

    const IdxT idx=vec.RowIdx->GetIdx();
    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
         sit!=mg.GetTriangVertexEnd(); ++sit){
        if (sit->Unknowns.Exist(idx))
            vec.Data[ sit->Unknowns(idx) ]= f(sit->GetCoord(),0.);
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
         sit!=mg.GetTriangEdgeEnd() && vec.RowIdx->NumUnknownsEdge(); ++sit){
        if (sit->Unknowns.Exist(idx))
            vec.Data[ sit->Unknowns(idx) ]= f(GetBaryCenter(*sit),0.);
    }
}

/// Set vectorial valued function on vertices and edges
void SetDOF(const DROPS::MultiGridCL& mg, DROPS::VecDescCL& vec, vec_fun f)
{
    using namespace DROPS;
    const IdxT idx=vec.RowIdx->GetIdx();
    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
         sit!=mg.GetTriangVertexEnd(); ++sit){
        if (sit->Unknowns.Exist(idx)){
            Point3DCL val= f(sit->GetCoord(), 0.);
            for (int i=0; i<3; ++i)
                vec.Data[ sit->Unknowns(idx)+i ]= val[i];
        }
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
         sit!=mg.GetTriangEdgeEnd() && vec.RowIdx->NumUnknownsEdge(); ++sit){
        if (sit->Unknowns.Exist(idx)){
            Point3DCL val= f(GetBaryCenter(*sit), 0.);
            for (int i=0; i<3; ++i)
                vec.Data[ sit->Unknowns(idx)+i ]= val[i];
        }
    }
}

/// Check if scalar values are read in correctly
double CheckDOF(const DROPS::MultiGridCL& mg, DROPS::VecDescCL& vec, scal_fun f)
{
    using namespace DROPS;
    double error=-1.;
    int counter=0;

    const IdxT idx=vec.RowIdx->GetIdx();

    // check data on vertices
    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
         sit!=mg.GetTriangVertexEnd(); ++sit){
        counter++;
        if (sit->Unknowns.Exist(idx)){
            const double diff=std::fabs(vec.Data[ sit->Unknowns(idx) ] - f(sit->GetCoord(),0.));

            if(diff > 0 )
            {
                std::cout << counter << ": Diff-Error : ("<< sit->GetGID() <<")"<<"Level:" << sit->GetLevel() <<std::endl;
                std::cout << sit->GetCoord() << " R: ["<<vec.Data[ sit->Unknowns(idx)] <<"]" << (double) (f(sit->GetCoord(),0.))<< std::endl;
            }

            if (diff>error)
                error= diff;
        }

    }

    // check data on edges
    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
         sit!=mg.GetTriangEdgeEnd() && vec.RowIdx->NumUnknownsEdge(); ++sit){
        if (sit->Unknowns.Exist(idx)){
            const double diff=std::fabs(vec.Data[ sit->Unknowns(idx) ] - f(GetBaryCenter(*sit),0.));
            if (diff>error)
                error= diff;
        }
    }

    return ProcCL::GlobalMax(error);
}

/// Check if vectorial values are read in correctly
bool CheckDOF(const DROPS::MultiGridCL& mg, DROPS::VecDescCL& vec, vec_fun f)
{
    using namespace DROPS;
    double error=-1.;

    const IdxT idx=vec.RowIdx->GetIdx();
    for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
         sit!=mg.GetTriangVertexEnd(); ++sit){
        if (sit->Unknowns.Exist(idx)){
            const Point3DCL val=MakePoint3D(vec.Data[ sit->Unknowns(idx)+0 ],
                                            vec.Data[ sit->Unknowns(idx)+1 ],
                                            vec.Data[ sit->Unknowns(idx)+2 ]);
            const double diff=(val-f(sit->GetCoord(),0.)).norm();
            if (diff>error)
                error= diff;
        }
    }
    for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
         sit!=mg.GetTriangEdgeEnd() && vec.RowIdx->NumUnknownsEdge(); ++sit){
        if (sit->Unknowns.Exist(idx)){
            const Point3DCL val=MakePoint3D(vec.Data[ sit->Unknowns(idx)+0 ],
                                            vec.Data[ sit->Unknowns(idx)+1 ],
                                            vec.Data[ sit->Unknowns(idx)+2 ]);
            const double diff=(val-f(GetBaryCenter(*sit),0.)).norm();
            if (diff>error)
                error= diff;
        }
    }
    return ProcCL::GlobalMax(error);
}


namespace DROPS{
/// \name Boundary of pressure and levelset
//@{
typedef BndDataCL<double> PressureBndCL;
BndCondT natural_cond[6]= {Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC};
PressureBndCL pr_bnd(6, natural_cond);
PressureBndCL lset_bnd(6, natural_cond);
//@}

/// \name Boundary of velocity
//@{
typedef BndDataCL<Point3DCL> VelocityBndCL;
BndCondT vel_cond[6]= {Dir0BC, Dir0BC, Dir0BC, DirBC, Dir0BC, Dir0BC};
const BndSegDataCL<Point3DCL>::bnd_val_fun bnd_fun_vel[6]={ZeroVel, ZeroVel, ZeroVel, ::Velocity, ZeroVel, ZeroVel};
VelocityBndCL vel_bnd(6, vel_cond, bnd_fun_vel);
//@}

/// \name FE-Types
//@{
typedef P1EvalCL<double,    const PressureBndCL, const VecDescCL> PrFE;
typedef P2EvalCL<double,    const PressureBndCL, const VecDescCL> LsetFE;
typedef P2EvalCL<Point3DCL, const VelocityBndCL, const VecDescCL> VelFE;
//@}
}

/****************************************************************************
* S E R I A L I Z A T I O N                                                 *
****************************************************************************/

// Write out the parallel multigrid and dof (if C.unknowns!=0) into files
void PerformSerialization(DROPS::LoadBalHandlerCL& lb, DROPS::ParMultiGridCL& pmg)
{
    using namespace DROPS;
    MultiGridCL& mg=lb.GetMG();

    if (DROPS::ProcCL::IamMaster())
        std::cout << line << " Refine the grid :\n"
                  << "  - "<<C.markall    <<" times regular\n"
                  << "  - "<<C.markdrop   <<" times around a droplet in the middle\n"
                  << "  - "<<C.markcorner <<" times in the origin corner" << std::endl;

    // Refine MG as given in the parameter-file
    for (int ref=0; ref<C.markall + C.markdrop + C.markcorner; ++ref){

        if (ref > 0 && ref < C.markall)
            MarkAll(mg);
        else if (C.markall <= ref && C.markall+C.markdrop > ref)
            MarkDrop(mg);
        else
            MarkCorner(mg);

        pmg.Refine();
        lb.DoMigration();
    }

    // Check multigrid, if everything is OK
    if (DROPS::ProcCL::IamMaster())
        std::cout << " Check the parallel multigrid ... " << flush;

    if (CheckParMultiGrid(pmg)){
        if (DROPS::ProcCL::IamMaster()) std::cout << " OK\n";
    }
    else {
        if (DROPS::ProcCL::IamMaster())
            std::cout << "  !!! error in MG !!!\n";
    }

    if (DROPS::ProcCL::IamMaster())
        cout << " Distribution of elements:\n";
    mg.SizeInfo(cout);

    if (DROPS::ProcCL::IamMaster())
        cout << line << " Write out multigrid\n" << line << flush;

    // Write out multigrid
    DROPS::ParMGSerializationCL writer(mg, C.ser_dir, DROPS::ProcCL::Master());
    writer.WriteMG();


    VTKOutCL *vtkwriter=0;
    if (C.vtk){
        vtkwriter = new VTKOutCL(mg, C.vtkName, 1, string(C.vtkDir + "/serialization") , C.vtkBinary);
        vtkwriter->PutGeom(0.0);
    }


    if (C.unknowns){
        // Create dof on vertices and edges
        IdxDescCL pr_idx  ( P1_FE);
        IdxDescCL lset_idx( P2_FE);
        IdxDescCL vel_idx ( vecP2_FE);

        pr_idx.CreateNumbering( mg.GetLastLevel(), mg, pr_bnd);
        lset_idx.CreateNumbering( mg.GetLastLevel(), mg, lset_bnd);
        vel_idx.CreateNumbering( mg.GetLastLevel(), mg, vel_bnd);

        VecDescCL pr_vec(&pr_idx);
        VecDescCL lset_vec(&lset_idx);
        VecDescCL vel_vec(&vel_idx);

        SetDOF(mg, pr_vec,   Pressure);
        SetDOF(mg, lset_vec, Levelset);
        SetDOF(mg, vel_vec,  Velocity);

        // write out dof
        if (DROPS::ProcCL::IamMaster())
            cout << " Write out dof into files:\n";

        writer.WriteDOF(&lset_vec, "Levelset");
        writer.WriteDOF(&pr_vec, "Pressure");
        writer.WriteDOF(&vel_vec, "Velocity");

        if (C.vtk){
            vtkwriter->PutScalar("Pressure", PrFE  (&pr_vec,   &pr_bnd,   &mg));
            vtkwriter->PutScalar("Levelset", LsetFE(&lset_vec, &lset_bnd, &mg));
            vtkwriter->PutVector("Velocity", VelFE (&vel_vec,  &vel_bnd,  &mg));
            vtkwriter->Commit();
        }
        pr_idx.DeleteNumbering( mg);
        lset_idx.DeleteNumbering( mg);
        vel_idx.DeleteNumbering( mg);
    }

     if (vtkwriter) delete vtkwriter;
}

void CheckSerialization(DROPS::LoadBalHandlerCL& lb, DROPS::ParMultiGridCL& pmg)
// Read a serialized multigrid, distribute the multigrid and check if everything
// is fine
{
    using namespace DROPS;
    MultiGridCL& mg= lb.GetMG();

    if (DROPS::ProcCL::IamMaster())
        cout << line << " Check if read multigrid is sane on master" << endl;

    if (DROPS::ProcCL::IamMaster()){
        ofstream serSanity("sanity.txt");
        serSanity << SanityMGOutCL(mg) << std::endl;
    }

    IdxDescCL *pr_idx=0, *lset_idx=0, *vel_idx=0;       // indices of read dof
    VecDescCL *pr_vec=0, *lset_vec=0, *vel_vec=0;       // vector describers of read dof

    // Read unknowns
    if (C.unknowns){
        pr_idx  = new IdxDescCL( P1_FE);
        lset_idx= new IdxDescCL( P2_FE);
        vel_idx = new IdxDescCL( vecP2_FE);

        if (ProcCL::IamMaster()){
            pr_idx->CreateNumbering( mg.GetLastLevel(), mg, pr_bnd);
            lset_idx->CreateNumbering( mg.GetLastLevel(), mg, lset_bnd);
            vel_idx->CreateNumbering( mg.GetLastLevel(), mg, vel_bnd);
        }

        pr_vec  = new VecDescCL(pr_idx);
        lset_vec= new VecDescCL(lset_idx);
        vel_vec = new VecDescCL(vel_idx);

        pmg.AttachTo(pr_vec,   &pr_bnd);
        pmg.AttachTo(lset_vec, &lset_bnd);
        pmg.AttachTo(vel_vec,  &vel_bnd);

        if (ProcCL::IamMaster()){
            ReadDOF(mg, pr_vec,   std::string(C.ser_dir+"Pressure"));
            ReadDOF(mg, lset_vec, std::string(C.ser_dir+"Levelset"));
            ReadDOF(mg, vel_vec,  std::string(C.ser_dir+"Velocity"));
        }

    }

    if (DROPS::ProcCL::IamMaster())
        cout << line << " Distribute MultiGrid ... "<< endl;

    lb.DoMigration();

    IdxDescCL *dist_pr_idx=0, *dist_lset_idx=0, *dist_vel_idx=0;       // indices of distributed dof
    VecDescCL *dist_pr_vec=0, *dist_lset_vec=0, *dist_vel_vec=0;       // vector describers of distributed dof

    if (C.unknowns){
        dist_pr_idx  = new IdxDescCL( P1_FE);
        dist_lset_idx= new IdxDescCL( P2_FE);
        dist_vel_idx = new IdxDescCL( vecP2_FE);

        dist_pr_idx->CreateNumbering( mg.GetLastLevel(), mg, pr_bnd);
        dist_lset_idx->CreateNumbering( mg.GetLastLevel(), mg, lset_bnd);
        dist_vel_idx->CreateNumbering( mg.GetLastLevel(), mg, vel_bnd);

        dist_pr_vec  = new VecDescCL(dist_pr_idx);
        dist_lset_vec= new VecDescCL(dist_lset_idx);
        dist_vel_vec = new VecDescCL(dist_vel_idx);

        pmg.HandleNewIdx(pr_idx, dist_pr_vec);
        pmg.HandleNewIdx(lset_idx, dist_lset_vec);
        pmg.HandleNewIdx(vel_idx, dist_vel_vec);
    }

    // Check multigrid, if everything is OK
    if (DROPS::ProcCL::IamMaster())
        std::cout << " Check the parallel multigrid ... " << flush;

    if (CheckParMultiGrid(pmg)){
        if (DROPS::ProcCL::IamMaster()) std::cout << " OK\n";
    }
    else {
        if (DROPS::ProcCL::IamMaster())
            std::cout << "  !!! error in MG !!!\n";
    }

    if (C.unknowns){
        // Check values of dof
        if (DROPS::ProcCL::IamMaster())
            std::cout << " Check values of dof ... " << std::endl;

        const double pr_error  = CheckDOF(mg, *dist_pr_vec,   Pressure);
        const double lset_error= CheckDOF(mg, *dist_lset_vec, Levelset);
        const double vel_error = CheckDOF(mg, *dist_vel_vec,  Velocity);

        if (DROPS::ProcCL::IamMaster())
            std::cout << "  - Greatest difference for pressure: "<< pr_error   <<'\n'
                      << "  - Greatest difference for levelset: "<< lset_error <<'\n'
                      << "  - Greatest difference for velocity: "<< vel_error  <<std::endl;
    }
    else{
        // Check if a refinement works ...
        if (DROPS::ProcCL::IamMaster())
            std::cout << " Check if a refinement and a load balancing step can be performed ... ";

        DROPS::MarkAll(mg);
        pmg.Refine();
        lb.DoMigration();

        if (DROPS::ProcCL::IamMaster())
            std::cout << "OK"<<std::endl;
    }
}

int main (int argc, char** argv)
{

    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
    try
    {
        if (argc<2){
            if (DROPS::ProcCL::IamMaster()){
              std::cerr << "usage: "<<argv[0]<<" <param-file> ";
            }
            throw DROPS::DROPSErrCL("No enough parameters are given!");
        }

        // Read configuration
        std::ifstream param( argv[1]);
        if (!param){
          std::cerr << "error while opening parameter file\n"; return 1;
        }
        param >> C;

        param.close();
        if (DROPS::ProcCL::IamMaster())
          std::cout << C << std::endl;

        // Init of the parallel structurs.
        DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

        if (DROPS::ProcCL::IamMaster() && C.mode==0)
            std::cout << line << " Create initial grid and distribution ... \n";
        if (DROPS::ProcCL::IamMaster() && C.mode==1)
            std::cout << line << " Read multigrid out of a file ... \n";

        // Create a unitcube at origin 0
        DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
        e1[0]=C.dim[0]; e2[1]=C.dim[1]; e3[2]= C.dim[2];

        DROPS::MGBuilderCL        * mgb     = 0;
        DROPS::BrickBuilderCL     * builder = 0;

        if (DROPS::ProcCL::IamMaster()){
            // if mode==0, then create a brick and write it out
            if (C.mode==0){
                mgb = new DROPS::BrickBuilderCL(C.orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
            }
            else{
                builder = new DROPS::BrickBuilderCL(C.orig, e1, e2, e3, C.basicref_x, C.basicref_y, C.basicref_z);
                mgb = new DROPS::FileBuilderCL(C.ser_dir, builder);
            }
        }
        else{
            mgb = new DROPS::EmptyBrickBuilderCL(C.orig, e1, e2, e3);
        }

        // Create multigrid and attach to pmg
        DROPS::MultiGridCL mg (*mgb) ;
        pmg.AttachTo(mg);
        DROPS::LoadBalHandlerCL lb(mg);

        // Do initial distribution of the elements

        if (C.mode==0){
            lb.DoInitDistribution(DROPS::ProcCL::Master());
            PerformSerialization(lb, pmg);
        }
        else{
            CheckSerialization(lb, pmg);
        }

        return 0;
    }
    catch (DROPS::DROPSErrCL err) { err.handle(); }

}


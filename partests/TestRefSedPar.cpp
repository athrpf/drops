//**************************************************************************
// File:    TestRefSedPar.cpp                                              *
// Content: Check refinement and load balancing for sedimenting droplet    *
//          without unknowns                                               *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   20. October 2008                                               *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestRefSedPar.cpp
/// \brief Check refinement and load balancing for sedimenting droplet without unknowns

// include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"
#include <ddd.h>

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include numeric computing!
#include "levelset/adaptriang.h"

 // include in- and output
#include "partests/params.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#ifdef __SUNPRO_CC
#  include <math.h>     // for pi
#endif


DROPS::ParamParBrickFlowCL C;

const char line[] ="------------------------------------------------------------";

void DisplayDetailedGeom(DROPS::MultiGridCL& mg)
{
    const DROPS::Uint level=mg.GetLastLevel();
    DROPS::Uint *numTetrasAllProc=0;
    DROPS::Uint *numFacesAllProc=0;
    DROPS::Uint *numDistFaceAllProc=0;
    if (DROPS::ProcCL::IamMaster()){
        numTetrasAllProc  = new DROPS::Uint[DROPS::ProcCL::Size()];
        numFacesAllProc   = new DROPS::Uint[DROPS::ProcCL::Size()];
        numDistFaceAllProc= new DROPS::Uint[DROPS::ProcCL::Size()];
    }
    // Gather information about distribution on master processor
    DROPS::Gather(mg.GetNumTriangTetra(level),      numTetrasAllProc,   DROPS::ProcCL::Master());
    DROPS::Gather(mg.GetNumTriangFace(level),       numFacesAllProc,    DROPS::ProcCL::Master());
    DROPS::Gather(mg.GetNumDistributedFaces(level), numDistFaceAllProc, DROPS::ProcCL::Master());

    // Display information
    if (DROPS::ProcCL::IamMaster()){
        double ratioTetra       =  (double)*std::max_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size())
                                  /(double)*std::min_element(numTetrasAllProc,   numTetrasAllProc+DROPS::ProcCL::Size());
        DROPS::Uint allTetra    =  std::accumulate(numTetrasAllProc, numTetrasAllProc+DROPS::ProcCL::Size(), 0),
                    allFace     =  std::accumulate(numFacesAllProc, numFacesAllProc+DROPS::ProcCL::Size(), 0),
                    allDistFace =  std::accumulate(numDistFaceAllProc, numDistFaceAllProc+DROPS::ProcCL::Size(), 0);
        double      *ratioDistFace=new double[DROPS::ProcCL::Size()];
        double maxRatio= *std::max_element(ratioDistFace, ratioDistFace+DROPS::ProcCL::Size());

        // global information
        std::cerr << "Detailed information about the parallel multigrid:\n"
                  << "#(master tetras on finest level):    "<<allTetra<<'\n'
                  << "#(all Faces on finest level):        "<<allFace<<'\n'
                  << "#(distributed Faces on fines level): "<<allDistFace<<'\n';
        std::cerr << "Ratio between max/min Tetra: "<<ratioTetra
                  <<" max Ratio DistFace/AllFace: "<<maxRatio<<std::endl;

        // local information for all processors
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            ratioDistFace[i]= ((double)numDistFaceAllProc[i]/(double)numFacesAllProc[i]*100.);

        std::cerr << std::setw(6)  <<  "Proc"
                  << std::setw(8)  << "#Tetra"
                  << std::setw(8)  << "#Faces"
                  << std::setw(12) << "#DistFaces"
                  << std::setw(12) << "%DistFaces"
                  << '\n';
        for (int i=0; i<DROPS::ProcCL::Size(); ++i)
            std::cerr << std::setw(6)  << i
                      << std::setw(8)  << numTetrasAllProc[i]
                      << std::setw(8)  << numFacesAllProc[i]
                      << std::setw(12) << numDistFaceAllProc[i]
                      << std::setw(12) << ratioDistFace[i] << std::endl;

        // free memory
        if (numTetrasAllProc)   delete[] numTetrasAllProc;
        if (numFacesAllProc)    delete[] numFacesAllProc;
        if (numDistFaceAllProc) delete[] numDistFaceAllProc;
        if (ratioDistFace)      delete[] ratioDistFace;
    }
}

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)
{ return DROPS::SVectorCL<3>(0.); }

namespace DROPS{
class MovingDropletCL
{
private:
    static double    SedVel_;
    static double    rad_;
    static Point3DCL dropletCenter_;

public:
    MovingDropletCL() {}

    /// \brief Init
    void Init(double vel, double rad, const Point3DCL& center){
        SedVel_=vel;
        rad_=rad;
        dropletCenter_=center;
    }

    /// \brief Distance of a point to the interface
    static double DistFunction( const DROPS::Point3DCL& p){
        const Point3DCL vectorToCenter   = dropletCenter_-p;
        return vectorToCenter.norm()-rad_;
    }

    /// \brief Move droplet up
    void Update() { dropletCenter_[1]+=SedVel_; }

    /// \brief Get center of the droplet
    Point3DCL GetCenter() const { return dropletCenter_; }

    /// \brief Get sedimentation velocity of the droplet
    double GetSedVel() const { return SedVel_; }
};

double    MovingDropletCL::SedVel_=double();
double    MovingDropletCL::rad_=double();
Point3DCL MovingDropletCL::dropletCenter_=Point3DCL();
}

namespace DROPS{

void Strategy( ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    MultiGridCL& mg= pmg.GetMG();
    ParTimerCL time;
    double durationRef, durationWrite;

    MovingDropletCL movingDroplet;
    movingDroplet.Init(C.Anstroem , C.Radius[0], C.Mitte);
    AdapTriangCL adapt( pmg, lb, C.ref_width, 0, C.ref_flevel);
    VTKOutCL vtkwriter(mg, C.vtkName, (C.vtk ? C.num_steps/C.vtk+1 : 0), 
                        std::string(C.vtkDir + "/" + C.vtkName), false);
    for (int i=1; i<=C.num_steps; ++i){
        if (ProcCL::IamMaster()){
            std::cerr << line << " Step "<<i<<std::endl;
            std::cerr << "Center of the droplet: "<<movingDroplet.GetCenter()<<std::endl;
            std::cerr << "Sedimentation velocity: "<<movingDroplet.GetSedVel()<<"("<<C.Anstroem<<")"<<std::endl;
        }
        time.Reset();
        adapt.MakeInitialTriang(MovingDropletCL::DistFunction);
        time.Stop(); durationRef=time.GetTime();
        
        time.Reset();
        if (C.vtk && i%C.vtk==0){
            vtkwriter.PutGeom(i*C.dt);
            vtkwriter.Commit();
        }
        time.Stop(); durationWrite=time.GetTime();
        if (ProcCL::IamMaster()){
            std::cerr << "Refinements for step "<<i<<" took: "<<durationRef<<" s.\n"
                      << "Writing geometry for setp "<<i<<" took "<<durationWrite<<"s."<<std::endl;
        }

        movingDroplet.Update();
    }

    //DisplayDetailedGeom(mg);
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
  try
  {
    if (argc!=2)
    {
        IF_MASTER
          std::cerr << "You have to specify one parameter:\n\t"
                    << argv[0] << " <param_file>" << std::endl;
        return 1;
    }
    std::ifstream param( argv[1]);
    if (!param)
    {
        IF_MASTER
          std::cerr << "error while opening parameter file\n";
        return 1;
    }
    param >> C;
    param.close();
    IF_MASTER
      std::cerr << C << std::endl;

    DROPS::ParMultiGridCL pmg(3);

    int nx, ny, nz;
    double dx, dy, dz;
    std::string mesh( C.meshfile), delim("x@");
    size_t idx;
    while ((idx= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
    if (!brick_info || dx!=dz)
    {
        std::cerr << "error while reading geometry information: " << mesh << "\n";
        return 1;
    }

    C.r_inlet= dx/2;
    DROPS::Point3DCL orig, px, py, pz;
    px[0]= dx; py[1]= dy; pz[2]= dz;


    DROPS::ParTimerCL time;
    DROPS::MGBuilderCL    *mgb     = 0;
    DROPS::BrickBuilderCL *builder = 0;
    if (DROPS::ProcCL::IamMaster()){
        if (C.IniCond!=3){
            mgb = new DROPS::BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
        }
        else{   // read geometry out of a file
            builder = new DROPS::BrickBuilderCL(orig, px, py, pz, nx, ny, nz);
            mgb = new DROPS::FileBuilderCL(C.ser_dir, builder);
        }
    }
    else{
        mgb = new DROPS::EmptyBrickBuilderCL(orig, px, py, pz);
    }

    DROPS::MultiGridCL mg(*mgb);
    pmg.AttachTo(mg);

    DROPS::LoadBalHandlerCL lb(mg);
    if (C.IniCond!=3)
        lb.DoInitDistribution(DROPS::ProcCL::Master());

    switch (C.refineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
    }

    Strategy( pmg, lb);    // do all the stuff

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


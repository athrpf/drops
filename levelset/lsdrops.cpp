//**************************************************************************
// File:    lsdrops.cpp                                                    *
// Content: test case for drop in stationary flow                          *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include <fstream>
#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "out/ensightOut.h"

DROPS::Point3DCL FixedVel( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL midpt( 0.5), q= p-midpt; 
    double d= q.norm(),
           c= d<0.25 ? d : (d<0.5 ? 0.5-d: 0);
    q[2]= q[0]; q[0]= c*q[1]; q[1]= -c*q[2]; q[2]= 0.; 
    return q;
}

double DistanceFct( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL midpt( 0.5);
    midpt[1]= 0.25;
    return (midpt-p).norm()-0.2;
}


namespace DROPS
{  // for strategy

template<class ProblemT>
void Strategy( ProblemT& prob, double dt, int num_steps, double SD, int num_reparam)
{
    MultiGridCL& mg= prob.GetMG();

    LevelsetP2CL<ProblemT> lset( mg, 0, 0.5, SD);
    IdxDescCL& lidx= lset.idx;
    IdxDescCL& vidx= prob.vel_idx;
    VecDescCL& vel=  prob.v;

    prob.CreateNumberingVel( mg.GetLastLevel(), &vidx);
    vel.SetIdx( &vidx);
    { // init vel
        IdxT pos= 0;
        for (MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(), end= mg.GetTriangVertexEnd();
            it!=end; ++it)
        {
            const Point3DCL val= FixedVel( it->GetCoord());
            vel.Data[pos++]= val[0]; 
            vel.Data[pos++]= val[1]; 
            vel.Data[pos++]= val[2]; 
        }
        for (MultiGridCL::TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(), end= mg.GetTriangEdgeEnd();
            it!=end; ++it)
        {
            const Point3DCL val= FixedVel( GetBaryCenter( *it));
            vel.Data[pos++]= val[0]; 
            vel.Data[pos++]= val[1]; 
            vel.Data[pos++]= val[2]; 
        }
    }

    lset.CreateNumbering( mg.GetLastLevel(), &lidx);
    lset.Phi.SetIdx( &lidx);
    lset.Init( DistanceFct);
    lset.SetupSystem( prob.GetVelSolution() );
    lset.SetTimeStep( dt);
    
    EnsightP2SolOutCL ensight( mg, &lidx);
    
    const char datgeo[]= "ensight/drop.geo", 
               datvec[]= "ensight/drop.vec",
               datscl[]= "ensight/drop.scl";
    ensight.CaseBegin( "drop.case", num_steps+1);
    ensight.DescribeGeom( "rotating vel field", datgeo);
    ensight.DescribeScalar( "Levelset", datscl, true); 
    ensight.DescribeVector( "Velocity", datvec); 
    ensight.putGeom( datgeo);
    ensight.putScalar( datscl, lset.GetSolution(), 0);
    ensight.putVector( datvec, prob.GetVelSolution());

    for (int i=1; i<=num_steps; ++i)
    {
        lset.DoStep();
        
        if (i%20==0)
            lset.Reparam( num_reparam, 0.01);
        // after half of the time, the velocity field is turned around
        if (i==num_steps/2)
        {
            prob.v.Data*= -1;
            lset.SetupSystem( prob.GetVelSolution() );
            lset.SetTimeStep( dt);
        }
        
        ensight.putScalar( datscl, lset.GetSolution(), i/10.);
    }

    ensight.CaseEnd();
}

class DummyStokesCoeffCL {};

} // end of namespace DROPS

int main( int argc, char **argv)
{
    double dt= 0.1, SD= 0.1;
    int num_steps= 100, num_reparam= 0;
    
    if (argc>1)
        dt= atof( argv[1]);
    if (argc>2)
        num_steps= atoi( argv[2]);
    if (argc>3)
        SD= atof( argv[3]);
    if (argc>4)
        num_reparam= atoi( argv[4]);
        
    std::cout << "dt = " << dt << ", SD = " << SD << ", num_reparam = " << num_reparam << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    typedef DROPS::StokesP2P1CL<DROPS::DummyStokesCoeffCL> StokesOnBrickCL;
    typedef StokesOnBrickCL                                                  MyStokesCL;
    
    int num;
    std::cout << "# Unterteilungen: "; std::cin >> num;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, num, num, num);
    const bool IsNeumann[6]= 
        {true, true, true, true, true, true};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &FixedVel, &FixedVel, &FixedVel, &FixedVel, &FixedVel, &FixedVel};
        
    MyStokesCL prob(brick, DROPS::DummyStokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
 
    Strategy( prob, dt, num_steps, SD, num_reparam);

    return 0;
}

//**************************************************************************
// File:    reparam.cpp                                                    *
// Content: test reparametrization of the levelset function                *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include <fstream>
#include "stokes/stokes.h"
#include "levelset/levelset.h"
#include "out/ensightOut.h"

DROPS::Point3DCL Null( const DROPS::Point3DCL&, double)
{
    return DROPS::Point3DCL(0.);
}

double SmoothedSign( double x, double alpha)
{
    return x/std::sqrt(x*x+alpha*alpha);
}

double DistFct( const DROPS::Point3DCL& p)
{
    DROPS::Point3DCL midpt( 0.5);
    midpt[1]= 0.25;
    return (midpt-p).norm()-0.2;
}

double Phi0( const DROPS::Point3DCL& p)
{
    double dist= DistFct( p);
    return SmoothedSign( dist, 0.1);
}

double Phi1( const DROPS::Point3DCL& p)
{
    return 6*DistFct( p);
}

double Phi2( const DROPS::Point3DCL& p)
{
    const double depth= 0.8-p[1],
                 alpha= 0.1;

    if (depth>0)
        return SmoothedSign( DistFct( p), alpha);
    else
        return SmoothedSign( depth, alpha);
}


namespace DROPS
{  // for strategy

template<class ProblemT>
void Strategy( ProblemT& prob, double dt, int num_steps, double diff, int bsp)
{
    MultiGridCL& mg= prob.GetMG();
    // Levelset: sigma= 0, theta= 0.5, SD= 0
    LevelsetP2CL lset( mg, 0, 0.5, 0, diff);

    IdxDescCL& lidx= lset.idx;
    IdxDescCL& vidx= prob.vel_idx;
    VecDescCL& vel=  prob.v;

    prob.CreateNumberingVel( mg.GetLastLevel(), &vidx);
    lset.CreateNumbering(    mg.GetLastLevel(), &lidx);
    vel.SetIdx( &vidx);
    lset.Phi.SetIdx( &lidx);
    switch (bsp)
    {
        case 0:  lset.Init( Phi0); break;
        case 1:  lset.Init( Phi1); break;
        case 2:  lset.Init( Phi2); break;
        default: lset.Init( DistFct);
    }
    lset.SetupSystem( prob.GetVelSolution() );

    EnsightP2SolOutCL sol( mg, &lidx);

    const char datgeo[]= "ensight/rep.geo",
               datscl[]= "ensight/rep.scl";
    sol.CaseBegin( "rep.case", num_steps+1);
    sol.DescribeGeom( "Cube", datgeo);
    sol.DescribeScalar( "Levelset", datscl, true);
    sol.putGeom( datgeo);
    sol.putScalar( datscl, lset.GetSolution(), 0.);

    TimerCL time;
    time.Start();
    lset.ReparamFastMarching();
    time.Stop();
    std::cerr << time.GetTime() << " sec for Fast Marching\n";
    sol.putScalar( datscl, lset.GetSolution(), dt/2);

    for (int i=1; i<=num_steps; ++i)
    {
        lset.Reparam( 1, dt);
        sol.putScalar( datscl, lset.GetSolution(), i*dt);
    }

    sol.CaseEnd();
}

class DummyStokesCoeffCL {};

} // end of namespace DROPS

int main( int argc, char **argv)
{
  try{
    double dt= 0.01, diff= 1e-4;
    int bsp= 0, num_steps= 100;

    if (argc>1)
        diff= std::atof( argv[1]);
    if (argc>2)
        dt= std::atof( argv[2]);
    if (argc>3)
        num_steps= std::atoi( argv[3]);

    std::cout << num_steps << " steps of lenght dt = " << dt
              << ", diff = " << diff << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    typedef DROPS::StokesP2P1CL<DROPS::DummyStokesCoeffCL> StokesOnBrickCL;
    typedef StokesOnBrickCL                                                  MyStokesCL;

    int num= 0;
    std::cout << "# Unterteilungen: "; std::cin >> num;
    std::cout << "Beispielnr.: "; std::cin >> bsp;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, num, num, num);
    const bool IsNeumann[6]=
        {true, true, true, true, true, true};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]=
        { &Null, &Null, &Null, &Null, &Null, &Null };

    MyStokesCL prob(brick, DROPS::DummyStokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));

    Strategy( prob, dt, num_steps, diff, bsp);

    return 0;
  } catch(DROPS::DROPSErrCL err) { err.handle(); }
}

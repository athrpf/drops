//**************************************************************************
// File:    reparam.cpp                                                    *
// Content: test reparametrization of the levelset function                *
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
void Strategy( ProblemT& prob, double dt, int num_steps, double SD, int bsp, int meth)
{
    MultiGridCL& mg= prob.GetMG();
    LevelsetP2CL lset( mg, 0, 1.0, SD);

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

    for (int i=1; i<=num_steps; ++i)
    {
/*
        if (meth)
            lset.Reparam2();
        else
*/
            lset.Reparam( 1, dt);
        sol.putScalar( datscl, lset.GetSolution(), i*dt);
    }
}

class DummyStokesCoeffCL {};

} // end of namespace DROPS

int main( int argc, char **argv)
{
    double dt= 0.1, SD= 0.1;
    int bsp, meth= 0, num_steps= 10;
    
    if (argc>1)
        dt= atof( argv[1]);
    if (argc>2)
        num_steps= atoi( argv[2]);
    if (argc>3)
        SD= atof( argv[3]);

    std::cout << num_steps << " steps of lenght dt = " << dt << ", SD = " << SD << std::endl;
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    typedef DROPS::StokesP2P1CL<DROPS::DummyStokesCoeffCL> StokesOnBrickCL;
    typedef StokesOnBrickCL                                                  MyStokesCL;
    
    int num;
    std::cout << "# Unterteilungen: "; std::cin >> num;
    std::cout << "Beispielnr.: "; std::cin >> bsp;
//    std::cout << "Reparam.verfahren (1=SaveIF, 0=Rep.Gl.) "; std::cin >> meth;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, num, num, num);
    const bool IsNeumann[6]= 
        {true, true, true, true, true, true};
    const DROPS::StokesVelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &FixedVel, &FixedVel, &FixedVel, &FixedVel, &FixedVel, &FixedVel};
        
    MyStokesCL prob(brick, DROPS::DummyStokesCoeffCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
 
    Strategy( prob, dt, num_steps, SD, bsp, meth);

    return 0;
}

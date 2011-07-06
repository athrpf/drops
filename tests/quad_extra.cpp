#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/discretize.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "num/quadrature.h"
#include <iomanip>

using namespace DROPS;

int degx, degy, degz;

double f(const Point3DCL& p, double)
{
    return  (degx==0 ? 1. : std::pow( p[0], degx))
           *(degy==0 ? 1. : std::pow( p[1], degy))
           *(degz==0 ? 1. : std::pow( p[2], degz));
}

double exactint[56] = { // with maple
  0.1666666666666666666666667, 0.04166666666666666666666667,
  0.01666666666666666666666667, 0.008333333333333333333333333,
  0.004761904761904761904761905, 0.002976190476190476190476190,
  0.04166666666666666666666667, 0.008333333333333333333333333,
  0.002777777777777777777777778, 0.001190476190476190476190476,
  0.0005952380952380952380952381, 0.01666666666666666666666667,
  0.002777777777777777777777778, 0.0007936507936507936507936508,
  0.0002976190476190476190476190, 0.008333333333333333333333333,
  0.001190476190476190476190476, 0.0002976190476190476190476190,
  0.004761904761904761904761905, 0.0005952380952380952380952381,
  0.002976190476190476190476190, 0.04166666666666666666666667,
  0.008333333333333333333333333, 0.002777777777777777777777778,
  0.001190476190476190476190476, 0.0005952380952380952380952381,
  0.008333333333333333333333333, 0.001388888888888888888888889,
  0.0003968253968253968253968254, 0.0001488095238095238095238095,
  0.002777777777777777777777778, 0.0003968253968253968253968254,
  0.00009920634920634920634920635, 0.001190476190476190476190476,
  0.0001488095238095238095238095, 0.0005952380952380952380952381,
  0.01666666666666666666666667, 0.002777777777777777777777778,
  0.0007936507936507936507936508, 0.0002976190476190476190476190,
  0.002777777777777777777777778, 0.0003968253968253968253968254,
  0.00009920634920634920634920635, 0.0007936507936507936507936508,
  0.00009920634920634920634920635, 0.0002976190476190476190476190,
  0.008333333333333333333333333, 0.001190476190476190476190476,
  0.0002976190476190476190476190, 0.001190476190476190476190476,
  0.0001488095238095238095238095, 0.0002976190476190476190476190,
  0.004761904761904761904761905, 0.0005952380952380952380952381,
  0.0005952380952380952380952381, 0.002976190476190476190476190
};

void TestExactness_extrapolation(int num_extrapolation)
{   
    DROPS::TetraBuilderCL tet( 0);
    DROPS::MultiGridCL mg( tet);
    TetraCL& s= *mg.GetAllTetraBegin();
//    s.DebugInfo( std::cout);
    std::cout.precision( 18);

    QuadDomainCL qdom;
    DROPS::ExtrapolationToZeroCL extra( num_extrapolation, DROPS::RombergSubdivisionCL());
    // DROPS::ExtrapolationToZeroCL extra( num_level, DROPS::HarmonicSubdivisionCL());
    LocalP2CL<> ones(-1.);
    make_ExtrapolatedQuad5Domain( qdom, ones, extra);
    //TetraPartitionCL partition;
    //partition.make_partition< SortedVertexPolicyCL,MergeCutPolicyCL>(2, ones);
    //make_CompositeQuad5Domain( qdom, partition);
    
    DROPS::GridFunctionCL<> integrand;
    
    Quad5CL<> q;
    size_t c= 0;
    for (degz= 0; degz <= 5; ++degz) {
        for (degy= 0; degy + degz <= 5; ++degy) {
            for (degx= 0; degx + degy + degz <= 5; ++degx) {
                //q.assign( s, f);
                resize_and_evaluate_on_vertexes (f, s, qdom, 0., integrand);
                std::cout << "degz: " << degz << "\tdegy: " << degy << "\tdegx: " << degx
                          << "\t\tI-Q_h: " << exactint[c++] - quad(integrand, 1., qdom, NegTetraC)//q.quad( 1.)
                          << "\tIntegral: " << quad(integrand, 1., qdom, NegTetraC) <<  "             ";//q.quad( 1.)
                /*for (size_t i= 0; i < q.size(); ++i)
                    std::cout << '\t' << q[i];*/
                std::cout << std::endl;
            }
        }
    }
}
inline double tetracut (const DROPS::Point3DCL& p)
{
    double s = 1.;
    if (p[2] == 0)
        s = 0.;
    // return p.norm() - 0.5;
    return s;
}

inline double tetracut_instat (const DROPS::Point3DCL& p, double)
{
    return tetracut( p);
}

int fakultaet (int n)
{
    int f =1;
    for (int it = 2; it<=n; ++it)
        f *=it;
    return f;
}

int binomi (int n, int k)
{
    return fakultaet(n)/(fakultaet(k)*fakultaet(n-k));
}
   

void TestExactness_extrapolation2D(int num_extrapolation)
{  
    double exactint_surf[21];     
    DROPS::TetraBuilderCL tet( 0);
    DROPS::MultiGridCL mg( tet);
    TetraCL& s= *mg.GetAllTetraBegin();
//    s.DebugInfo( std::cout);
    std::cout.precision( 18);

    DROPS::QuadDomain2DCL qdom;
    DROPS::ExtrapolationToZeroCL extra( num_extrapolation, DROPS::RombergSubdivisionCL());
    DROPS_FOR_TRIANG_TETRA( mg, 0, it) {
        DROPS::LocalP2CL<> ls_loc( *it,&tetracut_instat);
        make_ExtrapolatedQuad5Domain2D( qdom, ls_loc,*it, extra);
        DROPS::GridFunctionCL<> integrand;
        Quad5CL<> q;
        degz = 0;
        int it =0.;
        for (degy= 0; degy <= 5; ++degy) {
            for (degx= 0; degx + degy <= 5; ++degx) {
                exactint_surf[it] = 0.;
                for (int k = 0; k<= degy+1; ++k) {
                    exactint_surf[it] += std::pow(-1,k)*binomi(degy+1,k)/(degx+k+1);
                }    
                exactint_surf[it] *= 1./(degy+1);
                resize_and_evaluate_on_vertexes (f, s, qdom, 0., integrand);
                std::cout << "degz: " << degz << "\tdegy: " << degy << "\tdegx: " << degx
                          << "\t\tI-Q_h: " << exactint_surf[it] - quad_2D(integrand, qdom)//q.quad( 1.)
                          << "\tIntegral: " << quad_2D(integrand,  qdom) <<  "             ";//q.quad( 1.)
                /*for (size_t i= 0; i < q.size(); ++i)
                    std::cout << '\t' << q[i];*/
                std::cout << std::endl;
                ++it;
             }
        }
    }
}

int main ()
{
  try {
    TestExactness_extrapolation( 5);
    TestExactness_extrapolation2D( 5);
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
    return 0;
}
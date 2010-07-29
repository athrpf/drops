/// \file reparam.cpp
/// \brief test reparametrization of the levelset function
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

/*
 * This file is part of DROPS.
 *
 * DROPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DROPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with DROPS. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Copyright 2009 LNM/SC RWTH Aachen, Germany
*/

#include <fstream>
#include "levelset/levelset.h"
#include "levelset/adaptriang.h"
#include "levelset/params.h"
#include "levelset/surfacetension.h"
#include "out/ensightOut.h"
#include "out/vtkOut.h"
#include "misc/problem.h"

namespace DROPS{

/// \brief Class for determining distances to an ellipsoid
class EllipsoidCL
{
  private:
    static Point3DCL Mitte_;
    static Point3DCL Radius_;

  public:
    EllipsoidCL( const Point3DCL& Mitte, const Point3DCL& Radius)
        { Init( Mitte, Radius); }
    static void Init( const Point3DCL& Mitte, const Point3DCL& Radius)
        { Mitte_= Mitte;    Radius_= Radius; }
    static double DistanceFct( const Point3DCL& p)
    {
        Point3DCL d= p - Mitte_;
        const double avgRad= cbrt(Radius_[0]*Radius_[1]*Radius_[2]);
        d/= Radius_;
        return std::abs( avgRad)*d.norm() - avgRad;
    }
    static Point3DCL& GetCenter() { return Mitte_; }
    static Point3DCL& GetRadius() { return Radius_; }
};
Point3DCL EllipsoidCL::Mitte_;
Point3DCL EllipsoidCL::Radius_;

/// \brief Class for determining distances to a horizontal sliced domain
class HorizontalSlicesCL
{
  private:
    static int numSlices_;
    static double yBottom_;
    static double yTop_;

  public:
    HorizontalSlicesCL( int numSlices, double yBottom, double yTop)
        { Init( numSlices, yBottom, yTop); }
    static void Init( int numSlices, double yBottom, double yTop)
        { numSlices_=numSlices; yBottom_=yBottom; yTop_=yTop; }
    static double DistanceFct( const Point3DCL& p)
    {
        const double sizeOfSlice= (yTop_-yBottom_)/(numSlices_+1);
        double distance= std::numeric_limits<double>::max();
        for ( int i=1; i<=numSlices_; ++i){
            const double posOfSlice=yBottom_+i*sizeOfSlice;
            distance= std::min( distance, std::abs(p[1]-posOfSlice));
        }
        return distance;
    }
};
int HorizontalSlicesCL::numSlices_;
double HorizontalSlicesCL::yBottom_;
double HorizontalSlicesCL::yTop_;

/// \brief Signed distance function to a torus, negative in the torus.
/** The symmetry plane is \f$\{p_2=0\}\f$. The intersection with the
    torus is the annulus beween two circles centered in the origin and radii
    \f$R-r\f$, \f$R+r\f$. The intersections with planes containing the
    \f$e_2\f$-axis are composed of two circles with radius \f$r\f$, the
    centers of which are at distance \f$R\f$ from the origin.
*/
class TorusCL
{
  private:
    static double r_;
    static double R_;

public:
    TorusCL( double r, double R)
        { Init(r,R); }
    static void Init( double r, double R)
        { r_=r; R_=R; }
    static double DistanceFct( const Point3DCL& p)
    {
        return std::sqrt( p[2]*p[2] + std::pow( std::sqrt( p[0]*p[0] + p[1]*p[1]) - R_, 2) ) - r_;
    }
};
double TorusCL::r_=0.1;
double TorusCL::R_=0.3;

class MyParamCL: public ParamBrickCL,    public ParamReparamCL,
                 public ParamAdaptRefCL, public ParamExperimentalDataCL,
                 public ParamVTKCL
{
  protected:
    void RegisterParams() {}
  public:
    MyParamCL() { RegisterParams(); }
    MyParamCL( const string& filename) {
        std::ifstream file(filename.c_str());
        rp_.ReadParams( file);
        ParamBrickCL::rp_.ReadParams( file);
        ParamReparamCL::rp_.ReadParams( file);
        ParamAdaptRefCL::rp_.ReadParams( file);
        ParamExperimentalDataCL::rp_.ReadParams( file);
        ParamVTKCL::rp_.ReadParams( file);
    }
} C;

double sigmaf (const Point3DCL&, double) { return 0.; }
Point3DCL gsigma (const Point3DCL&, double) { return Point3DCL(); }

/// \brief Disturb a given level set function, i.e., scale vector with 100.
void Disturb( VectorCL& phi){
    phi *=100.;
}

/** \brief Check result of the re-parametrization
    Here, the following errors besides gradient information are gathered of the reparametrizated level set function:
    \f{eqnarray*}{
     e_1 &=& \frac{1}{|S|} \sum_{v\in S}|\phi^{\mbox{ex}}(v) - \phi^{\mbox{comp}}(v)|\\
     e_1^{\mbox{rel}}&=& \frac{1}{|S|} \sum_{v\in S} \frac{|\phi^{\mbox{ex}}(v) - \phi^{\mbox{comp}}(v)|}{|\phi^{\mbox{ex}}(v)|}\\
     e_2&=& \sqrt{\frac{1}{|S|} \sum_{v\in S}|\phi^{\mbox{ex}}(v) - \phi^{\mbox{comp}}(v)|^2}\\
     e_2^{\mbox{rel}}&=& \sqrt{\frac{1}{|S|} \sum_{v\in S} \frac{|\phi^{\mbox{ex}}(v) - \phi^{\mbox{comp}}(v)|^2}{|\phi^{\mbox{ex}}(v)|^2}}\\
     e_{\infty}&=& \max_{v\in S}|\phi^{\mbox{ex}}(v) - \phi^{\mbox{comp}}(v)|\\
     e_{\infty}^{\mbox{rel}}&=&  \max_{v\in S} \frac{|\phi^{\mbox{ex}}(v) - \phi^{\mbox{comp}}(v)|}{|\phi^{\mbox{ex}}(v)|},
     }
    where S denotes the set of vertices not lying at the interface, phi^comp
    the computed, and phi^ex the given level set function.
*/
void CheckReparametrization( const LevelsetP2CL& lset, const VectorCL& phiEx)
{
    double maxGradPhi, minGradPhi;
    double e1=0., e1Rel=0., e2=0., e2Rel=0., eSup=0., eSupRel=0.;  // error on off-site vertices

    lset.GetMaxMinGradPhi( maxGradPhi, minGradPhi);

    int   num_sign[3];
    IdxT  Numb[10];
    int   sign[10];
    const Uint idx= lset.Phi.RowIdx->GetIdx();
    const RefRuleCL RegRef= GetRefRule( RegRefRuleC);
    std::vector<bool> frontier( phiEx.size(), false);

    DROPS_FOR_TRIANG_CONST_TETRA( lset.GetMG(), -1, it){
        for ( int v=0; v<10; ++v){ // collect data on all DoF
            if (v<4)
                Numb[v]= it->GetVertex(v)->Unknowns(idx);
            else
                Numb[v]= it->GetEdge(v-4)->Unknowns(idx);
            sign[v]= std::abs( lset.Phi.Data[Numb[v]]<1e-8 ? 0 : ( lset.Phi.Data[Numb[v]]>0 ? 1 : -1));
            if ( lset.Phi.Data[Numb[v]]<1e-8)
                frontier[Numb[v]]= true;
        }

        // Check if child tetrahedron is intersected
        for ( Uint ch=0; ch<MaxChildrenC; ++ch){
            const ChildDataCL data= GetChildData( RegRef.Children[ch]);
            num_sign[0]= num_sign[1]= num_sign[2]= 0;
            for ( Uint vert= 0; vert<NumVertsC; ++vert)
                ++num_sign[ sign[data.Vertices[vert]] + 1];
            const bool intersec= (num_sign[0]*num_sign[2]!=0); // Vorzeichenwechsel
            // child is intersected by interface
            if ( intersec){
                for ( Uint vert=0; vert<4; ++vert){
                    frontier[ Numb[data.Vertices[vert]]]=true;
                }
            }
        }
    }

    // Compute errors
    DROPS_FOR_TRIANG_CONST_VERTEX( lset.GetMG(), -1, it){
#ifdef _PAR
        if ( !it->IsExclusive()) continue;
#endif
        IdxT dof= it->Unknowns( idx);
        if (!frontier[dof]){
            const double phiRep = lset.Phi.Data[dof], phiExact= phiEx[ dof];
            const double diff   = std::abs( phiRep-phiExact);
            const double relDiff= diff/std::abs(phiExact);
            e1     += diff;
            e2     += diff*diff;
            e1Rel  += relDiff;
            e2Rel  += relDiff*relDiff;
            eSup    = std::max( eSup, diff);
            eSupRel = std::max( eSupRel, relDiff);
        }
    }
    DROPS_FOR_TRIANG_CONST_EDGE( lset.GetMG(), -1, it){
#ifdef _PAR
        if ( !it->IsExclusive()) continue;
#endif
        IdxT dof= it->Unknowns( idx);
        if (!frontier[dof]){
            const double phiRep = lset.Phi.Data[dof], phiExact= phiEx[ dof];
            const double diff   = std::abs( phiRep-phiExact);
            const double relDiff= diff/std::abs(phiExact);
            e1     += diff;
            e2     += diff*diff;
            e1Rel  += relDiff;
            e2Rel  += relDiff*relDiff;
            eSup    = std::max( eSup, diff);
            eSupRel = std::max( eSupRel, relDiff);
        }
    }

    Uint numOffsite=0;
    for ( size_t i=0; i< frontier.size(); ++i) if ( !frontier[i]) ++numOffsite;
#ifdef _PAR
    numOffsite= ProcCL::GlobalSum( numOffsite, ProcCL::Master());
    e1        = ProcCL::GlobalSum( e1, ProcCL::Master());
    e1Rel     = ProcCL::GlobalSum( e1Rel, ProcCL::Master());
    e2        = ProcCL::GlobalSum( e2, ProcCL::Master());
    e2Rel     = ProcCL::GlobalSum( e2Rel, ProcCL::Master());
    eSup      = ProcCL::GlobalMax( eSup, ProcCL::Master());
    eSupRel   = ProcCL::GlobalMax( eSupRel, ProcCL::Master());
#endif
    e1    /= numOffsite;
    e1Rel /= numOffsite;
    e2     = std::sqrt( e2/numOffsite);
    e2Rel  = std::sqrt( e2Rel/numOffsite);

#ifdef _PAR
    if ( lset.Phi.RowIdx->GetEx().IsAcc( lset.Phi.Data))
        std::cout << " Vector is accumulated!" << std::endl;
    else
        std::cout << " Vector is not accumulated!" << std::endl;
#endif

    std::cout << "\n----------------------------------------------\n"
              << "Difference of re-parametrized level set function to given level set function:\n"
              << "e1 " << e1 << ", e1 (rel) " << e1Rel
              << ", e2 " << e2 << ", e2 (rel) " << e2Rel
              << ", eSup " << eSup << ", eSup (rel) " << eSupRel
              << "\nmax and min grad: " << maxGradPhi << ' ' << minGradPhi << std::endl;
}

void Strategy( DROPS::AdapTriangCL& adap, DROPS::BndDataCL<>& lsbnd)
{
    SurfaceTensionCL sf( sigmaf, gsigma);   // dummy class

    LevelsetP2CL lset( adap.GetMG(), lsbnd, sf);

    // writer for vtk-format
    VTKOutCL vtkwriter(adap.GetMG(), "DROPS data", (C.vtk_VTKOut ? 3 : 0),
                std::string(C.vtk_VTKDir + "/" + C.vtk_VTKName), C.vtk_Binary);
    vtkwriter.Register( make_VTKScalar( lset.GetSolution(), "level-set") );

    // Create numbering and assign given distance function
    lset.CreateNumbering( adap.GetMG().GetLastLevel(), &lset.idx);
    lset.Phi.SetIdx( &lset.idx);

    // Write out information:
    size_t numLsetUnk= lset.Phi.Data.size();
#ifdef _PAR
    numLsetUnk= ProcCL::GlobalSum( numLsetUnk);
#endif
    std::cout << numLsetUnk << " (accumulated) levelset unknowns.\n\n";

    switch ( C.rpm_Freq){
        case -1 : 
            std::cout << "Taking torus of radi (" << C.exp_RadDrop[0] << ',' << C.exp_RadDrop[1] << ") as level set function\n" << std::endl;
            lset.Init( TorusCL::DistanceFct); 
            break;
        case  0 : 
            std::cout << "Taking ellipsoid at " << C.exp_PosDrop << " and radi " << C.exp_RadDrop << " as level set function\n" << std::endl;
            lset.Init( EllipsoidCL::DistanceFct); 
            break;
        default:  
            std::cout << "Taking " << C.rpm_Freq << " horizontal sclices as level set function" << std::endl;
            lset.Init( HorizontalSlicesCL::DistanceFct);
    }

    if (C.vtk_VTKOut){
        vtkwriter.Write(0.0, true);
    }

    // Make a copy of exact level set function
    VectorCL phiEx( lset.Phi.Data);

    // Disturb level set function
    Disturb( lset.Phi.Data);

    // Perform re-parametrization
    std::auto_ptr<ReparamCL> reparam= ReparamFactoryCL::GetReparam( adap.GetMG(), lset.Phi, C.rpm_Method, /*periodic*/ false, &lset.GetBndData());
    reparam->Perform();

//    FastMarchCL fmm( adap.GetMG(), lset.Phi);
//    fmm.Reparam( true, 1);

    if (C.vtk_VTKOut){
        vtkwriter.Write(1.0, true);
    }

    // Check result
    CheckReparametrization( lset, phiEx);

    // Write difference as output
    VectorCL phiDiff( lset.Phi.Data-phiEx);
    if (C.vtk_VTKOut){
        std::swap( lset.Phi.Data, phiDiff);
        vtkwriter.Write(0.0, true);
        std::swap( lset.Phi.Data, phiDiff);
    }

}

} // end of namespace DROPS

int main( int argc, char **argv)
{
#ifdef _PAR
    DROPS::ProcInitCL procinit(&argc, &argv);
#endif
    try {
#ifdef _PAR
        DROPS::ParMultiGridInitCL pmginit;
#endif

        using DROPS::C;
        std::ifstream param;
        if (argc != 2) {
            std::cout << "Using default parameter file: reparam.param\n";
            param.open("reparam.param");
        }
        else{
            std::cout << "Opening file " << argv[1] << std::endl;
            param.open(argv[1]);
        }
        if (!param) {
            std::cerr << "error while opening parameter file\n";
            return 1;
        }
        param >> C;
        param.close();
        std::cout << DROPS::C << std::endl;

        DROPS::MultiGridCL* mg= 0;
        DROPS::BrickBuilderCL *mgb = 0;
        DROPS::Point3DCL a,b,c;
        a[0]= C.brk_dim[0];
        b[1]= C.brk_dim[1];
        c[2]= C.brk_dim[2];
        IF_MASTER
            mgb = new DROPS::BrickBuilderCL( C.brk_orig, a, b, c, C.brk_BasicRefX, C.brk_BasicRefY, C.brk_BasicRefZ);
        IF_NOT_MASTER
            mgb = new DROPS::EmptyBrickBuilderCL(C.brk_orig, a, b, c);

        mg= new DROPS::MultiGridCL( *mgb);
        delete mgb;

        DROPS::AdapTriangCL adap( *mg, C.ref_Width, C.ref_CoarsestLevel, C.ref_FinestLevel, -1);

        DROPS::EllipsoidCL::Init( C.exp_PosDrop, C.exp_RadDrop);
        DROPS::HorizontalSlicesCL::Init( C.rpm_Freq, C.brk_orig[1], C.brk_orig[1]+C.brk_dim[1] );
        DROPS::TorusCL::Init( C.exp_RadDrop[0], C.exp_RadDrop[1]);
        typedef double (*distance_fct)(const DROPS::Point3DCL&);
        distance_fct distance= C.rpm_Freq>0 ? DROPS::HorizontalSlicesCL::DistanceFct : DROPS::EllipsoidCL::DistanceFct;
        adap.MakeInitialTriang( distance);

        const DROPS::BndCondT bcls[6]= { DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC, DROPS::NoBC };
        const DROPS::LsetBndDataCL::bnd_val_fun bfunls[6]= { 0,0,0,0,0,0};
        DROPS::LsetBndDataCL lsbnd( 6, bcls, bfunls);

        DROPS::Strategy( adap, lsbnd);


    } catch (DROPS::DROPSErrCL err) {
        err.handle();
    }
}

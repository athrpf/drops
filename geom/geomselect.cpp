/// \file geomselect.h
/// \brief offers build/create routines for some standard domains
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Yuanjun Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#include "geom/geomselect.h"
#include <sstream>

namespace DROPS {

void BuildDomain( MultiGridCL* &mgp, const std::string& meshfile_name, int GeomType, const std::string& deserialization_file, double& r_inlet, std::vector<BndCondT>& BC)
{
#ifdef _PAR
    ParMultiGridCL::InstancePtr();
#endif
    if (GeomType == 0) {
        std::ifstream meshfile( meshfile_name.c_str());
        if (!meshfile)
            throw DROPSErrCL ("error while opening mesh file\n");

        ReadMeshBuilderCL *mgb= 0;       // builder of the multigrid

        // read geometry information from a file and create the multigrid
        IF_MASTER
            mgb = new ReadMeshBuilderCL( meshfile );
        IF_NOT_MASTER
            mgb = new EmptyReadMeshBuilderCL( meshfile );
        // Create the multigrid
        if (deserialization_file == "none")
            mgp= new MultiGridCL( *mgb);
        else {
            FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new MultiGridCL( filebuilder);
        }
        mgb->GetBC( BC);
        delete mgb;
    }
    if (GeomType == 1) {
        int nx, ny, nz;
        double dx, dy, dz;
        std::string mesh( meshfile_name), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz;
        if (!brick_info)
            throw DROPSErrCL("error while reading geometry information: " + mesh);
        r_inlet= dx/2;
        Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;

        BrickBuilderCL *mgb = 0;
        IF_MASTER
            mgb = new BrickBuilderCL( orig, px, py, pz, nx, ny, nz);
        IF_NOT_MASTER
            mgb = new EmptyBrickBuilderCL(orig, px, py, pz);

        if (deserialization_file == "none")
            mgp= new MultiGridCL( *mgb);
        else {
            FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new MultiGridCL( filebuilder);
        }
        delete mgb;
    }
    if (GeomType == 2) {
        int nx, ny, nz;
        double dx, dy, dz;
        int cnx, cny, cnz;
        Uint cdx, cdy, cdz;
        std::string mesh( meshfile_name), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz >> cdx >> cdy >> cdz >> cnx >> cny >> cnz;
        if (!brick_info)
            throw DROPSErrCL("error while reading geometry information: " + mesh);
        r_inlet= dx/2;
        Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;

        CavityBuilderCL *mgb = 0;
        SArrayCL<Uint, 3> corg, cav;
        corg[0]= cdx; corg[1]= cdy; corg[2]= cdz;
        cav[0]= cnx; cav[1]= cny; cav[2]= cnz;
        IF_MASTER
            mgb = new CavityBuilderCL( orig, px, py, pz, nx, ny, nz, corg, cav);
        IF_NOT_MASTER
            mgb = new EmptyCavityBuilderCL(orig, px, py, pz, 1, corg, cav);

        if (deserialization_file == "none")
            mgp= new MultiGridCL( *mgb);
        else {
            FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new MultiGridCL( filebuilder);
        }
        delete mgb;
    }
#ifndef _PAR
    if (GeomType == 3) { //LBuilder
        ///TODO: SegFault in the Builder
        int nx, ny, nz;
        double dx, dy, dz;
        int bx,by;
        std::string mesh( meshfile_name), delim("x@");
        size_t idx;
        while ((idx= mesh.find_first_of( delim)) != std::string::npos )
            mesh[idx]= ' ';
        std::istringstream brick_info( mesh);
        brick_info >> dx >> dy >> dz >> nx >> ny >> nz >> bx >> by ;
        if (!brick_info)
            throw DROPSErrCL("error while reading geometry information: " + mesh);
        r_inlet= dx/2;
        Point3DCL orig, px, py, pz;
        px[0]= dx; py[1]= dy; pz[2]= dz;

        LBuilderCL *mgb = 0;
        mgb = new LBuilderCL( orig, px, py, pz, nx, ny, nz, bx, by);

        if (deserialization_file == "none")
            mgp= new MultiGridCL( *mgb);
        else {
            FileBuilderCL filebuilder( deserialization_file, mgb);
            mgp= new MultiGridCL( filebuilder);
        }
        delete mgb;
    }
    if (GeomType == 4) { //BBuilder
            int nx, ny, nz;
            double dx, dy, dz;
            int bx,by,bz;
            std::string mesh( meshfile_name), delim("x@");
            size_t idx;
            while ((idx= mesh.find_first_of( delim)) != std::string::npos )
                mesh[idx]= ' ';
            std::istringstream brick_info( mesh);
            brick_info >> dx >> dy >> dz >> nx >> ny >> nz >> bx >> by >> bz ;
            if (!brick_info)
                throw DROPSErrCL("error while reading geometry information: " + mesh);
            r_inlet= dx/2;
            Point3DCL orig, px, py, pz;
            px[0]= dx; py[1]= dy; pz[2]= dz;

            BBuilderCL *mgb = 0;
            mgb = new BBuilderCL( orig, px, py, pz, nx, ny, nz, bx, by, bz);

            if (deserialization_file == "none")
                mgp= new MultiGridCL( *mgb);
            else {
                FileBuilderCL filebuilder( deserialization_file, mgb);
                mgp= new MultiGridCL( filebuilder);
            }
            delete mgb;
        }
#endif
}

/// \brief helper functions that determine some of the boundaries of a poisson problem
//@{
inline double Grad2X(const DROPS::Point3DCL& p, double t)
  { return (std::exp(t)*std::exp(p[0]+p[1]+p[2])); }
inline double Grad2Y(const DROPS::Point3DCL& p, double t)
  { return (std::exp(t)*std::exp(p[0]+p[1]+p[2])); }
inline double Grad2Z(const DROPS::Point3DCL& p, double t)
  { return (std::exp(t)*std::exp(p[0]+p[1]+p[2])); }

inline double Grad3X(const DROPS::Point3DCL&, double)
  { return (4.0); }
inline double Grad3Y(const DROPS::Point3DCL&, double)
  { return (0.0); }
inline double Grad3Z(const DROPS::Point3DCL&, double)
  { return (0.0); }
//}@

/// \brief boundary description of a neumann problem
class NeuVal2CL
{
  public:
    template<int seg>
    static double neu_val(const DROPS::Point3DCL& p, double t)
    {
      switch (seg)
      {
        case 0: return (-Grad2X( p, t));
        case 1: return Grad2X( p, t);
        case 2: return (-Grad2Y( p, t));
        case 3: return Grad2Y( p, t);
        case 4: return (-Grad2Z( p, t));
        case 5: return Grad2Z( p, t);
        default:
        {
          std::cout <<"error: neu_val";
          return 1;
        }
      }
    }
};

/// \brief boundary description of a neumann problem
class NeuVal3CL
{
  public:
    template<int seg>
    static double neu_val(const DROPS::Point3DCL& p, double t)
    {
      switch (seg)
      {
        case 0: return (-Grad3X( p, t));
        case 1: return Grad3X( p, t);
        case 2: return (-Grad3Y( p, t));
        case 3: return Grad3Y( p, t);
        case 4: return (-Grad3Z( p, t));
        case 5: return Grad3Z( p, t);
        default:
        {
          std::cout <<"error: neu_val";
          return 1;
        }
      }
    }
};

/// \brief boundary description of a neumann problem
class NeuVal4CL
{
  public:
    // boundary functions (neumann, dirichlet type)
    // used for BndSegCL-object of a UnitCube
    template< int seg>
    static double neu_val(const DROPS::Point3DCL& p, double= 0.0) { return -64.0*p[0]*p[1]*(1.0-p[0])*(1.0-p[1]); }
};


void BuildPoissonBoundaryData( MultiGridCL* &mgp, PoissonBndDataCL* &bnddata,
        instat_scalar_fun_ptr fun, int GeomType, int bnd_type, std::vector<BndCondT>& BC)
{
    if (GeomType == 0) {
        const BoundaryCL& bnd= mgp->GetBnd();
        const BndIdxT num_bnd= bnd.GetNumBndSeg();

        BndCondT* bc = new BndCondT[num_bnd];
        DROPS::BndDataCL<>::bnd_val_fun* bnd_fun = new DROPS::BndDataCL<>::bnd_val_fun[num_bnd];
        for (BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= BC[ i])==DirBC ? fun : &Zero;
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
        }
        bnddata = new PoissonBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
    }
    if (GeomType == 1) { //Boundary condition for BrickBuilder
        BndCondT bc[6]= { DirBC, DirBC, DirBC, DirBC, DirBC, DirBC };
        DROPS::BndDataCL<>::bnd_val_fun bfun[6]=
            { fun, fun, fun, fun, fun, fun};
        switch ( bnd_type) {
            case 1 : break; //hom. Dirichlet boundary conditions
            case 2 : {
                bc[0] = bc[1] = bc[2] = bc[3] = bc[4] = bc[5] = NatBC;
                bfun[0] = &DROPS::NeuVal2CL::neu_val<0>;
                bfun[1] = &DROPS::NeuVal2CL::neu_val<1>;
                bfun[2] = &DROPS::NeuVal2CL::neu_val<2>;
                bfun[3] = &DROPS::NeuVal2CL::neu_val<3>;
                bfun[4] = &DROPS::NeuVal2CL::neu_val<4>;
                bfun[5] = &DROPS::NeuVal2CL::neu_val<5>;
            } break;
            case 3 : {
                bc[0] = bc[1] = bc[2] = bc[3] = bc[4] = bc[5] = NatBC;
                bfun[0] = &DROPS::NeuVal3CL::neu_val<0>;
                bfun[1] = &DROPS::NeuVal3CL::neu_val<1>;
                bfun[2] = &DROPS::NeuVal3CL::neu_val<2>;
                bfun[3] = &DROPS::NeuVal3CL::neu_val<3>;
                bfun[4] = &DROPS::NeuVal3CL::neu_val<4>;
                bfun[5] = &DROPS::NeuVal3CL::neu_val<5>;
            } break;
            case 4 : {
                bc[1] = bc[3] = bc[4] = bc[5] = Nat0BC;
                bc[0] = Dir0BC;
                bc[2] = NatBC;
                bfun[0] = &Zero;
                bfun[1] = &Zero;
                bfun[2] = fun;
                bfun[3] = &Zero;
                bfun[4] = &Zero;
                bfun[5] = &Zero;
            } break;
            case 5 : { //drops_neumann.cpp
                bc[0] = bc[1] = bc[2] = bc[3] = bc[4] = Dir0BC;
                bc[5] = NatBC;
                bfun[0] = &Zero;
                bfun[1] = &Zero;
                bfun[2] = &Zero;
                bfun[3] = &Zero;
                bfun[4] = &Zero;
                bfun[5] = &DROPS::NeuVal4CL::neu_val<5>;
            } break;
        }
        bnddata = new PoissonBndDataCL(6, bc, bfun);
    }
    if (GeomType == 2) {
        BndCondT bc[12]= { Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC };
        DROPS::BndDataCL<>::bnd_val_fun bfun[12]=
            { &Zero, &Zero, &Zero, &Zero, &Zero, &Zero, &Zero, &Zero, &Zero, &Zero, &Zero, &Zero };
        bnddata = new PoissonBndDataCL(12, bc, bfun);
    }
    if (GeomType == 3) {
        BndCondT bc[14]= { DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC};
        DROPS::BndDataCL<>::bnd_val_fun bfun[14]=
            { fun, fun, fun, fun, fun, fun, fun, &DROPS::Zero, fun, fun, fun, fun, fun, fun};
        bnddata = new PoissonBndDataCL(14, bc, bfun);
    }
    if (GeomType == 4) { //for BBuilder
        BndCondT bc[24]= { DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC,
                           DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC};
        DROPS::PoissonBndDataCL::bnd_val_fun bfun[24]=
            { fun, fun, fun, fun, fun, fun, fun, fun,
              fun, fun, fun, fun, fun, fun, fun, fun,
              fun, fun, fun, fun, fun, fun, fun, fun };
        bnddata = new PoissonBndDataCL(24, bc, bfun);
    }
}


void BuildStokesBoundaryData( MultiGridCL* &mgp, StokesBndDataCL* &bnddata,
        instat_vector_fun_ptr inflow, int GeomType, int bnd_type, std::vector<BndCondT>& BC)
{
    if (GeomType == 0) {
        const BoundaryCL& bnd= mgp->GetBnd();
        const BndIdxT num_bnd= bnd.GetNumBndSeg();

        BndCondT* bc = new BndCondT[num_bnd];
        StokesVelBndDataCL::bnd_val_fun* bnd_fun = new StokesVelBndDataCL::bnd_val_fun[num_bnd];
        for (BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= BC[ i])==DirBC ? inflow : &ZeroVel;
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
        }
        bnddata = new StokesBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
    }
    if (GeomType == 1) {
        BndCondT bc[6]= { Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC };
        StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[6]=
            { &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel };
        switch (bnd_type)
        {
            case 1 : // hom. Dirichlet bnd-data
              break;
            case 2 : // inhom. Dirichlet bnd-data for in-/outflow
            {
                bc[2]= bc[3]= DirBC;
                bfun[2]= bfun[3]= inflow;
            } break;
            case 3 : // tube/canal
            {
                bc[3]= DirBC;
                bc[2]= NatBC; //Rohr
                //bc[2]= bc[4]= bc[5]= NatBC;          //Kanal
                bfun[2]= &ZeroVel;
                //bfun[2]=bfun[4]=bfun[5]= &ZeroVel;   //Kanal
                bfun[3]= inflow;
            } break;
            case 4 : // channel (eindhoven)
            {
                bc[0]= DirBC;
                bc[1]= NatBC; //Rohr
                bfun[1]= &ZeroVel;
                bfun[0]= inflow;
            } break;
            case 5 : // predefined velocity over all boundaries
            {
                bc[0] = bc[1] = bc[2] = bc[3] = bc[4] = bc[5] = DirBC;
                bfun[0] = bfun[1] = bfun[2] = bfun[3] = bfun[4] = bfun[5] = inflow;
            } break;
            case 6 : //driven cavity
            {
                bc[5] = DirBC;
                bfun[5] = inflow;
            }break;
            default: throw DROPSErrCL("Unknown boundary data type");
        }
        bnddata = new StokesBndDataCL(6, bc, bfun);
    }
    if (GeomType == 2) {
        BndCondT bc[12]= { Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC };
        StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[12]=
            { &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel, &ZeroVel };
        bnddata = new StokesBndDataCL(12, bc, bfun);
    }
}

void BuildLsetBoundaryData( MultiGridCL* &mgp, LsetBndDataCL* &bnddata,
       __UNUSED__ instat_scalar_fun_ptr inflow, int GeomType, __UNUSED__ int bnd_type, __UNUSED__ std::vector<BndCondT>& BC)
{
    if (GeomType == 0) {
        const BoundaryCL& bnd= mgp->GetBnd();
        const BndIdxT num_bnd= bnd.GetNumBndSeg();

        BndCondT* bc = new BndCondT[num_bnd];
        LsetBndDataCL::bnd_val_fun* bnd_fun = new BndDataCL<>::bnd_val_fun[num_bnd];
        for (BndIdxT i=0; i<num_bnd; ++i)
        {
            bc[i]= NoBC;
            bnd_fun[i]= 0;
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
        }
        bnddata = new LsetBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
    }
    if (GeomType == 1) {
        BndCondT bc[6]= { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
        LsetBndDataCL::bnd_val_fun bfun[6]= {0,0,0,0,0,0};
        bnddata = new LsetBndDataCL(6, bc, bfun);
    }
    if (GeomType == 2) {
        BndCondT bc[12]= { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC, NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
        LsetBndDataCL::bnd_val_fun bfun[12]= {0,0,0,0,0,0,0,0,0,0,0,0};
        bnddata = new LsetBndDataCL(12, bc, bfun);
    }
}

/// \brief Create geometry of a Mzelle or a brick
void CreateGeom (MultiGridCL* &mgp, StokesBndDataCL* &bnddata, LsetBndDataCL* &lsetbnddata,
                 instat_vector_fun_ptr inflow, instat_scalar_fun_ptr lsetinflow,
                 const std::string& meshfile_name,
                 int GeomType, int bnd_type,
                 const std::string& deserialization_file, double& r_inlet)
{
    std::vector<BndCondT> BC;
    BuildDomain( mgp, meshfile_name, GeomType, deserialization_file, r_inlet, BC);
    BuildStokesBoundaryData( mgp, bnddata, inflow, GeomType, bnd_type, BC);
    BuildLsetBoundaryData( mgp, lsetbnddata, lsetinflow, GeomType, bnd_type, BC);
    std::cout << "Generated MG of " << mgp->GetLastLevel() << " levels." << std::endl;
}

void CreateGeomPoisson (MultiGridCL* &mgp, PoissonBndDataCL* &bnddata,
                 instat_scalar_fun_ptr inflow,
                 const std::string& meshfile_name,
                 int GeomType, int bnd_type,
                 const std::string& deserialization_file, double& r_inlet)
{
    std::vector<BndCondT> BC;
    BuildDomain( mgp, meshfile_name, GeomType, deserialization_file, r_inlet, BC);
    BuildPoissonBoundaryData( mgp, bnddata, inflow, GeomType, bnd_type, BC);
    std::cout << "Generated MG of " << mgp->GetLastLevel() << " levels." << std::endl;
}

} //end of namespace drops

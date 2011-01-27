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
#include "misc/bndmap.h"

namespace DROPS {


template<typename T>
void readBoundary (T& bndcond, const std::string& bnd_type, int size){

	size_t idx;
	std::string bndtype( bnd_type), delim("!");

	while ((idx= bndtype.find_first_of( delim)) != std::string::npos )
		bndtype[idx]= ' ';

	std::istringstream bnd( bndtype);

	for(int i=0; i<size;++i){
		bnd >> bndcond[i];
	}

}

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

void BuildPoissonBoundaryData( MultiGridCL* &mgp, PoissonBndDataCL* &bnddata,
        int GeomType, const std::string& bnd_type, const std::string& bnd_funcs, std::vector<BndCondT>& BC)
{
	instat_scalar_fun_ptr Zero = InScaMap::getInstance().find("Zero")->second;
    if (GeomType == 0) {

    	string bfunc[1];
    	readBoundary(bfunc, bnd_funcs, 1);
    	instat_scalar_fun_ptr fun= InScaMap::getInstance().find(bfunc[0])->second;
    	const BoundaryCL& bnd= mgp->GetBnd();
        const BndIdxT num_bnd= bnd.GetNumBndSeg();

        BndCondT* bc = new BndCondT[num_bnd];
        DROPS::BndDataCL<>::bnd_val_fun* bnd_fun = new DROPS::BndDataCL<>::bnd_val_fun[num_bnd];
        for (BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= BC[ i])==DirBC ? fun : InScaMap::getInstance().find("Zero")->second;
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
        }
        bnddata = new PoissonBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
    }
    if (GeomType == 1) { //Boundary condition for BrickBuilder

    	BndCondT bc[6];
        DROPS::BndDataCL<>::bnd_val_fun bfun[6];

    	int bcc[6];
    	string bfuncs[6];

    	readBoundary(bcc, bnd_type, 6);
    	readBoundary(bfuncs, bnd_funcs, 6);

    	for(int i=0; i<6;++i){
    	   	bc[i] = BndCondT(bcc[i]);
    	   	bfun[i] = InScaMap::getInstance().find(bfuncs[i])->second;
     	}

        bnddata = new PoissonBndDataCL(6, bc, bfun);
    }
    if (GeomType == 2) {
        BndCondT bc[12]= { Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC };
        DROPS::BndDataCL<>::bnd_val_fun bfun[12]=
            { Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero, Zero };
        bnddata = new PoissonBndDataCL(12, bc, bfun);
    }
    if (GeomType == 3) {

    	string bfunc[1];
    	readBoundary(bfunc, bnd_funcs, 1);
    	instat_scalar_fun_ptr fun= InScaMap::getInstance().find(bfunc[0])->second;

        BndCondT bc[14]= { DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC, DirBC};
        DROPS::BndDataCL<>::bnd_val_fun bfun[14]=
            { fun, fun, fun, fun, fun, fun, fun, Zero, fun, fun, fun, fun, fun, fun};
        bnddata = new PoissonBndDataCL(14, bc, bfun);
    }
    if (GeomType == 4) { //for BBuilder

     	string bfunc[1];
    	readBoundary(bfunc, bnd_funcs, 1);

    	instat_scalar_fun_ptr fun= InScaMap::getInstance().find(bfunc[0])->second;

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
        int GeomType, const std::string& bnd_funcs, const std::string& bnd_type, std::vector<BndCondT>& BC)
{
	StokesVelBndDataCL::bnd_val_fun ZeroVel = InVecMap::getInstance().find("ZeroVel")->second;
    if (GeomType == 0) {

    	string bfunc[1];
    	readBoundary(bfunc, bnd_funcs, 1);

    	StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun;
    	bfun=InVecMap::getInstance().find(bfunc[0])->second;

        const BoundaryCL& bnd= mgp->GetBnd();
        const BndIdxT num_bnd= bnd.GetNumBndSeg();

        BndCondT* bc = new BndCondT[num_bnd];
        StokesVelBndDataCL::bnd_val_fun* bnd_fun = new StokesVelBndDataCL::bnd_val_fun[num_bnd];

        for (BndIdxT i=0; i<num_bnd; ++i)
        {
            bnd_fun[i]= (bc[i]= BC[ i])==DirBC ? bfun : ZeroVel;
            std::cout << "Bnd " << i << ": "; BndCondInfo( bc[i], std::cout);
        }
        bnddata = new StokesBndDataCL(num_bnd, bc, bnd_fun);
        delete[] bc;
        delete[] bnd_fun;
    }
    if (GeomType == 1) {

    	DROPS::BndCondT bc[6];
    	StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[6];

  	    int bcc[6];
   	    std::string bfunc[6];

   	    readBoundary(bcc, bnd_type, 6);
   	    readBoundary(bfunc, bnd_funcs, 6);
   		for(int i=0; i<6;++i){
   			bc[i] = DROPS::BndCondT(bcc[i]);
   			bfun[i]=InVecMap::getInstance().find(bfunc[i])->second;
   		}

        bnddata = new StokesBndDataCL(6, bc, bfun);
    }
    if (GeomType == 2) {
        BndCondT bc[12]= { Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC, Dir0BC };
        StokesBndDataCL::VelBndDataCL::bnd_val_fun bfun[12]=
            { ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel, ZeroVel };
        bnddata = new StokesBndDataCL(12, bc, bfun);
    }
}

void BuildLsetBoundaryData( MultiGridCL* &mgp, LsetBndDataCL* &bnddata,
       __UNUSED__ instat_scalar_fun_ptr inflow, int GeomType, __UNUSED__ const std::string& bnd_type, __UNUSED__ std::vector<BndCondT>& BC)
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
                  instat_scalar_fun_ptr lsetinflow, const std::string& meshfile_name,
                  int GeomType,const std::string& bnd_funcs, const std::string& bnd_type,
                  const std::string& deserialization_file, double& r_inlet)
{
    std::vector<BndCondT> BC;
    BuildDomain( mgp, meshfile_name, GeomType, deserialization_file, r_inlet, BC);
    BuildStokesBoundaryData( mgp, bnddata, GeomType, bnd_funcs, bnd_type, BC);
    BuildLsetBoundaryData( mgp, lsetbnddata, lsetinflow, GeomType, bnd_type, BC);
    std::cout << "Generated MG of " << mgp->GetLastLevel() << " levels." << std::endl;
}

void CreateGeomPoisson (MultiGridCL* &mgp, PoissonBndDataCL* &bnddata,
                          const std::string& meshfile_name, int GeomType,
                          const std::string& bnd_type, const std::string& bnd_funcs,
                          const std::string& deserialization_file, double& r_inlet)
{
    std::vector<BndCondT> BC;
    BuildDomain( mgp, meshfile_name, GeomType, deserialization_file, r_inlet, BC);
    BuildPoissonBoundaryData( mgp, bnddata, GeomType, bnd_type, bnd_funcs, BC);
    std::cout << "Generated MG of " << mgp->GetLastLevel() << " levels." << std::endl;
}

} //end of namespace drops

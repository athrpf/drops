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

void BuildDomain( MultiGridCL* &mgp, const std::string& meshfile_name, int GeomType, const std::string& deserialization_file, double& r_inlet)
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

} //end of namespace drops

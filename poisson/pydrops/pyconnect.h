/// \file   pyconnect.h
/// \brief  interpret python data for Drops' use
/// \author LNM RWTH Aachen: Liang Zhang
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

#ifndef PYCONNECT_ME_H
#define PYCONNECT_ME_H

#include "misc/container.h"
#include "geom/multigrid.h"
#include "pdefunction.h"

#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <math.h>

#include <boost/shared_ptr.hpp>

/// holds the Python input matrices and Python output parameters
class PythonConnectCL
{
public:
  typedef boost::shared_ptr<PythonConnectCL> Ptr;
  int Nx_, Ny_, Nz_, Nt_, Nxy_, Nyz_, Nxz_, Nxyz_; // N=number of points
  double dx_, dy_, dz_, dt_;
  bool   adjoint_;
  std::ofstream* outfile;

  //const DropsFunction *presol_, *DelPsi_;
  double* C3D_,                               // output matrices: temp solution (Nxyz x nt),
    *MaxIter_;                                // max. iterations of solver (1 x 1)
  //helper maps for barycenters
  //static DROPS::MultiGridCL* MG_;
  DROPS::MultiGridCL* MG_;

  //static FACE_MAP face_map_;
  //static TETRA_MAP tetra_map_;
  FACE_MAP face_map_;
  TETRA_MAP tetra_map_;


  int GetNum( const DROPS::Point3DCL& p, double t, int seg) const
  {
    if (seg == 0 || seg == 1){//yoz
      return (rd(p[2]/dz_)*Ny_ + rd(p[1]/dy_) + rd(t/dt_)*Ny_*Nz_);
    }
    if (seg == 2 || seg == 3){//xoz
      return (rd(p[2]/dz_)*Nx_ + rd(p[0]/dx_) + rd(t/dt_)*Nx_*Nz_);
    }
    if (seg == 4 || seg == 5){//xoy
      return (rd(p[1]/dy_)*Nx_ + rd(p[0]/dx_) + rd(t/dt_)*Nx_*Ny_);
    }
  }
  //
  int GetNum( const DROPS::Point3DCL& p, double t=0.) const
  {
    return (rd(p[2]/dz_)*Nxy_ + rd(p[1]/dy_)*Nx_ + rd(p[0]/dx_) + rd(t/dt_)*Nxyz_);
  }

  bool GetNum(const DROPS::Point3DCL& p, double t, int& ix, int& iy, int& iz, int& it) const {
    ix = rd(p[0]/dx_);
    iy = rd(p[1]/dy_);
    iz = rd(p[2]/dz_);
    it = rd(t/dt_);
    }
 public:
  PythonConnectCL()
    {
      Nx_=Ny_=Nz_=Nt_=Nxy_=Nyz_=Nxz_=Nxyz_=-1;
      dx_=dy_=dz_=0.0;
      face_map_.clear();
      tetra_map_.clear();
    }
  void SetMG( DROPS::MultiGridCL* MG) { MG_= MG; }
  //
  void ClearMaps() {
    face_map_.clear();
    tetra_map_.clear();
  }
  //
  void setFaceMap()
  {
    *outfile<<"SETTING FACE MAP" <<std::endl;
    DROPS::Uint lvl= MG_->GetLastLevel();
    for(DROPS::MultiGridCL::TriangFaceIteratorCL
	  fit=MG_->GetTriangFaceBegin(lvl), fend=MG_->GetTriangFaceEnd(lvl);
        fit != fend;++fit) {
      DROPS::FaceCL& face = *fit;

      DROPS::Point3DCL bc = DROPS::GetBaryCenter(face);
      d_pair pr= std::make_pair(rnd(bc[2]), rnd(bc[1]));
      cmp_key key= std::make_pair(rnd(bc[0]), pr);

      face_map_[key]= &face;
    }
    *outfile<<"FACE MAP SET" <<std::endl;
  }
  void DumpFaceMap()
  {
    FILE* f=fopen("face.map","w");
    *outfile<<"DUMPING FACE MAP\n"<<std::endl;
    for (face_it p= face_map_.begin(); p!= face_map_.end(); p++) {

      DROPS::Point3DCL bc = DROPS::Point3DCL(0.);
      bc[0]=p->first.first;
      bc[1]=p->first.second.second;
      bc[2]=p->first.second.first;

      d_pair pr= std::make_pair(bc[2], bc[1]);
      cmp_key key= std::make_pair(bc[0], pr);

      DROPS::FaceCL* face = face_map_[key];

      fprintf(f,"bc: %f %f %f \n",bc[0],bc[1],bc[2]);

      fprintf(f,"face:\n");
      DROPS::Point3DCL p = face->GetVertex(0)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = face->GetVertex(1)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = face->GetVertex(2)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
    }
    fprintf(f," %d",(int)face_map_.size());
    fclose(f);
    *outfile<<"END DUMP FACE MAP"<<std::endl;
  }
  //
  void setTetraMap()
  {
    *outfile<<"SETTING TETRA MAP"<<std::endl;
    DROPS::Uint lvl= MG_->GetLastLevel();
    for(DROPS::MultiGridCL::TriangTetraIteratorCL
	  tit=MG_->GetTriangTetraBegin(lvl), tend=MG_->GetTriangTetraEnd(lvl);
        tit != tend; ++tit) {
      DROPS::TetraCL& tetra = *tit;

      DROPS::Point3DCL bc = DROPS::GetBaryCenter(tetra);
      d_pair pr= std::make_pair(rnd(bc[2]), rnd(bc[1]));
      cmp_key key= std::make_pair(rnd(bc[0]), pr);

      tetra_map_[key]= &tetra;
    }
    *outfile<<"TETRA MAP SET"<<std::endl;
  }
  void DumpTetraMap()
  {
    FILE* f=fopen("tetra.map","w");
    *outfile<<"DUMPING TETRA MAP"<<std::endl;
    for (tetra_it p= tetra_map_.begin(); p!= tetra_map_.end(); p++) {

      DROPS::Point3DCL bc = DROPS::Point3DCL(0.);
      bc[0]=p->first.first;
      bc[1]=p->first.second.second;
      bc[2]=p->first.second.first;

      d_pair pr= std::make_pair(bc[2], bc[1]);
      cmp_key key= std::make_pair(bc[0], pr);

      DROPS::TetraCL* tetra = tetra_map_[key];
      fprintf(f,"key: %f %f %f \n",bc[0],bc[1],bc[2]);

      fprintf(f,"tetra:\n");
      DROPS::Point3DCL p = tetra->GetVertex(0)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = tetra->GetVertex(1)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = tetra->GetVertex(2)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
      p = tetra->GetVertex(3)->GetCoord();
      fprintf(f," %f %f %f \n",p[0],p[1],p[2]);
    }
    fclose(f);
    *outfile<<"END DUMP TETRA MAP"<<std::endl;
  }

  DropsFunction::Ptr GetDiffusion;
  DropsFunction::Ptr GetInitial;
  DropsFunction::Ptr GetSource;
  DropsFunction::Ptr GetInterfaceValue; // never barycentric
  DropsFunction::Ptr GetInflow; // never barycentric
  DropsFunction::Ptr GetDelPsi;
  DropsFunction::Ptr GetPresol;

  template<class P1EvalT>
    void SetSol3D( const P1EvalT& sol, double t)  //Instationary problem
    {
      double *out;
      if (adjoint_) {
	  const int num = (Nt_-rd(t/dt_))*Nxyz_;    //flip solution back
	  out = C3D_+num;
      } else {
	  const int num= (rd(t/dt_)-1)*Nxyz_;         // omit initial time step in output
	  out= C3D_+num;                              //don't like this way
      }

      DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
	  {
        out[GetNum( sit->GetCoord())]= sol.val( *sit);
	  }

    }
  template<class P1EvalT>
    void SetSol3D( const P1EvalT& sol)          //Stationary problem
    {
        double *out= C3D_;

        DROPS_FOR_TRIANG_CONST_VERTEX( sol.GetMG(), sol.GetLevel(), sit)
        {
            out[GetNum( sit->GetCoord())]= sol.val( *sit);
        }
    }

  //Check the input matrices
  void Init( std::ofstream* of, const DROPS::ParamCL& P, PdeFunction::ConstPtr C0, PdeFunction::ConstPtr B_in, PdeFunction::ConstPtr F, PdeFunction::ConstPtr Dw, PdeFunction::ConstPtr B_Inter, double* c_sol)
  {
    outfile=of;
    int refinesteps_;
    double lx_, ly_, lz_;
    int nx_, ny_, nz_;
    refinesteps_= P.get<int>("DomainCond.RefineSteps");

    std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
    size_t idx_;
    while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx_]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> lx_ >> ly_ >> lz_ >> nx_ >> ny_ >> nz_;
    Nx_ = nx_ * pow (2, refinesteps_)+1;
    Ny_ = ny_ * pow (2, refinesteps_)+1;
    Nz_ = nz_ * pow (2, refinesteps_)+1;
    Nt_ = P.get<int>("Time.NumSteps")+1;
    Nyz_=Ny_*Nz_; Nxy_=Nx_*Ny_; Nxz_=Nx_*Nz_;
    Nxyz_= Nxy_*Nz_;
    Nt_    = P.get<int>("Time.NumSteps")+1;
    dx_= lx_/(Nx_-1); dy_= ly_/(Ny_-1); dz_= lz_/(Nz_-1);
    dt_    = P.get<double>("Time.StepSize");

    // Save the matrix input arguments.
    GridFunction::Ptr vg(new VolumeGridFunction(dx_, dy_, dz_, dt_, &tetra_map_));
    GridFunction::Ptr sg_inlet(new SurfaceGridFunction(dx_, dy_, dz_, dt_, &face_map_, 0));
    GridFunction::Ptr sg_interface(new SurfaceGridFunction(dx_, dy_, dz_, dt_, &face_map_, 3));

    GetInitial = DropsFunction::Ptr(new DropsFunction(C0, vg, 4));
    GetInflow = DropsFunction::Ptr(new DropsFunction(B_in, sg_inlet, 3));
    GetSource = DropsFunction::Ptr(new DropsFunction(F, vg, 4));
    GetDiffusion = DropsFunction::Ptr(new DropsFunction(Dw, vg, 4));
    GetInterfaceValue = DropsFunction::Ptr(new DropsFunction(B_Inter, sg_interface, 3));

    std::string adstr ("IA1Adjoint");
    std::string IAProbstr = P.get<std::string>("PoissonCoeff.IAProb");
    adjoint_ = (adstr.compare(IAProbstr)==0);
    // Set the output pointer to the output arguments.
    C3D_ = c_sol;
  }

  void Init( std::ofstream* of, const DROPS::ParamCL& P, PdeFunction::ConstPtr B_in, PdeFunction::ConstPtr B_Inter, PdeFunction::ConstPtr F,
                    PdeFunction::ConstPtr presol, PdeFunction::ConstPtr DelPsi,PdeFunction::ConstPtr Dw,  double* c_sol)
  {
    outfile=of;
    int refinesteps_;
    double lx_, ly_, lz_;
    int nx_, ny_, nz_;
    refinesteps_= P.get<int>("DomainCond.RefineSteps");

    std::string mesh( P.get<std::string>("DomainCond.MeshFile")), delim("x@");
    size_t idx_;
    while ((idx_= mesh.find_first_of( delim)) != std::string::npos )
        mesh[idx_]= ' ';
    std::istringstream brick_info( mesh);
    brick_info >> lx_ >> ly_ >> lz_ >> nx_ >> ny_ >> nz_;
    Nx_ = nx_ * pow (2, refinesteps_)+1;
    Ny_ = ny_ * pow (2, refinesteps_)+1;
    Nz_ = nz_ * pow (2, refinesteps_)+1;
    Nyz_=Ny_*Nz_; Nxy_=Nx_*Ny_; Nxz_=Nx_*Nz_;
    Nxyz_= Nxy_*Nz_;
    dx_= lx_/(Nx_-1); dy_= ly_/(Ny_-1); dz_= lz_/(Nz_-1);
    dt_     = P.get<double>("Time.StepSize");

    GridFunction::Ptr vg(new VolumeGridFunction(dx_, dy_, dz_, dt_, &tetra_map_));
    GridFunction::Ptr sg_inlet(new SurfaceGridFunction(dx_, dy_, dz_, dt_, &face_map_, 0));
    GridFunction::Ptr sg_interface(new SurfaceGridFunction(dx_, dy_, dz_, dt_, &face_map_, 3));

    GetInflow = DropsFunction::Ptr(new DropsFunction(B_in, sg_inlet, 3));
    GetInterfaceValue = DropsFunction::Ptr(new DropsFunction(B_Inter, sg_interface, 3));
    GetSource = DropsFunction::Ptr(new DropsFunction(F, vg, 4));
    GetPresol = DropsFunction::Ptr(new DropsFunction(presol, vg, 4));
    GetDelPsi = DropsFunction::Ptr(new DropsFunction(DelPsi, vg, 4));
    GetDiffusion = DropsFunction::Ptr(new DropsFunction(Dw, vg, 4));

    // Set the output pointer to the output arguments.
    C3D_ = c_sol;
  }
};


class PyDropsErr : public DROPS::DROPSErrCL {
public:
 PyDropsErr(const PythonConnectCL::Ptr PyC_, DROPS::ParamCL* P_, int it_, std::string msg_) : PyC(PyC_), P(P_), it(it_), msg(msg_) {}
  const PythonConnectCL::Ptr PyC;
  DROPS::ParamCL* P;
  int it;
  std::string msg;

  /// Write all information corresponding to this problem to a file and return the filename
  std::string write_err_to_file() const
  {
    std::string filename(tmpnam(NULL));
    std::ofstream outfile(filename.c_str());
    std::string ll("--------------------------------");
    outfile << "----------------------------\nOutput file generated from a PyDropsErr Exception\n----------------------------\n";
    outfile << "1. Param file\n";
    outfile << (*P) <<std::endl;
    outfile << ll << "\nThe solver failed in time step " << it << " of "<< PyC->Nt_-1 << std::endl;
    /*
    std::stringstream Fstream, C0stream, Binstream, Binterfacestream, Dwstream;
    for (int ix=0; ix<PyC->Nx_; ++ix)
      for (int iy=0; iy<PyC->Ny_; ++iy)
	for (int iz=0; iz<PyC->Nz_; ++iz) {
	  Fstream << "F(" << ix << ","<<iy<<","<<iz<<") = " << (*PyC->F_)(ix,iy,iz,it) << std::endl;
	  C0stream << "C0(" << ix << ","<<iy<<","<<iz<<") = " << (*PyC->C0_)(ix,iy,iz,it) << std::endl;
	  Binstream << "Bin(" << ix << ","<<iy<<","<<iz<<") = " << (*PyC->B_in_)(ix,iy,iz,it) << std::endl;
	  Binterfacestream << "Binter(" << ix << ","<<iy<<","<<iz<<") = " << (*PyC->B_Inter_)(ix,iy,iz,it) << std::endl;
	  Dwstream << "Dw(" << ix << ","<<iy<<","<<iz<<") = " << (*PyC->Dw_)(ix,iy,iz,it) << std::endl;
	}
    outfile << ll << "\n2. Source term\n" << Fstream.str();
    outfile << ll << "\n3. Initial Value\n" << C0stream.str();
    outfile << ll << "\n4. Inlet BC\n" << Binstream.str();
    outfile << ll << "\n5. Interface BC\n" << Binterfacestream.str();
    outfile << ll << "\n6. Diffusion Coefficient\n" << Dwstream.str();
    */
    outfile.close();
    return filename;
  }
private:
  PyDropsErr();
};

#endif

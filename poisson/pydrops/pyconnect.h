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
#include <string>
#include <sstream>
#include <math.h>

inline int rd( double d) { return static_cast<int>( d+0.5); }                   // rounding

class PythonConnectCL
{ // holds the Python input matrices and Python output parameters
  typedef std::pair<double, double> d_pair;
  typedef std::pair<double, d_pair> cmp_key;
  //
  typedef std::map<cmp_key, DROPS::FaceCL*>  FACE_MAP;
  typedef std::map<cmp_key, DROPS::TetraCL*> TETRA_MAP;

 private:
  int Nx_, Ny_, Nz_, Nxy_, Nyz_, Nxz_, Nxyz_; // N=number of points
  double dx_, dy_, dz_, dt_;
  double D_mol_;

  const PdeFunction *C0_, *B_in_, *B_Inter_, *F_,  // initial+boundary+rhs function,
    *Dw_;                                     // wavy induced diffusion parameter as a function,
  PdeFunction *f1_, *f2_;
  double* C3D_,                               // output matrices: temp solution (Nxyz x nt),
    *MaxIter_;                                // max. iterations of solver (1 x 1)
  //helper maps for barycenters
  static DROPS::MultiGridCL* MG_;

  static FACE_MAP face_map_;
  static TETRA_MAP tetra_map_;

  typedef FACE_MAP::const_iterator fi;
  typedef TETRA_MAP::const_iterator ti;


  static double rnd(double d) {// rounding four digits after comma
    int i_d = (int)(d*10000);
    double d_d = (double)i_d/10000;
    return d_d;
  }
  //
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
      Nx_=Ny_=Nz_=Nxy_=Nyz_=Nxz_=Nxyz_=-1;
      dx_=dy_=dz_=0.0;

      face_map_.clear();
      tetra_map_.clear();
    }
  void SetMG( DROPS::MultiGridCL* MG) { MG_= MG; }
  //
  static void ClearMaps() {
    face_map_.clear();
    tetra_map_.clear();
  }
  //
  static void setFaceMap()
  {
    std::cout<<"SETTING FACE MAP" <<std::endl;
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
    std::cout<<"FACE MAP SET" <<std::endl;
  }
  static void DumpFaceMap()
  {
    FILE* f=fopen("face.map","w");
    std::cout<<"DUMPING FACE MAP\n"<<std::endl;
    for (fi p= face_map_.begin(); p!= face_map_.end(); p++) {

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
    std::cout<<"END DUMP FACE MAP"<<std::endl;
  }
  //
  static void setTetraMap()
  {
    std::cout<<"SETTING TETRA MAP"<<std::endl;
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
    std::cout<<"TETRA MAP SET"<<std::endl;
  }
  static void DumpTetraMap()
  {
    FILE* f=fopen("tetra.map","w");
    std::cout<<"DUMPING TETRA MAP"<<std::endl;
    for (ti p= tetra_map_.begin(); p!= tetra_map_.end(); p++) {

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
    std::cout<<"END DUMP TETRA MAP"<<std::endl;
  }
  //

  double GetProductF1( const DROPS::Point3DCL& p, double t)
  {
    int ix, iy, iz, it;
    GetNum(p,t,ix,iy,iz,it);
    return (*f1_)(ix,iy,iz,it);
  };

  double GetProductF2( const DROPS::Point3DCL& p, double t)
  {
    int ix, iy, iz, it;
    GetNum(p,t,ix,iy,iz,it);
    return (*f2_)(ix,iy,iz,it);
  };

   double GetInitial( const DROPS::Point3DCL& p, double t)
  {
    t=0.;
    int ix, iy, iz, it;
    GetNum(p,t,ix,iy,iz,it);
    return (*C0_)(ix,iy,iz,it);
  };
  //boundary functions
  //x=0;
   double GetInflow( const DROPS::Point3DCL& p, double t)
  {
    int ix, iy, iz, it;
    GetNum(p,t,ix,iy,iz,it);
    return (*B_in_)(ix,iy,iz,it);
  };
  //y=0: if neumann condition is active
   double GetInterfaceFlux( const DROPS::Point3DCL& p, double t)
  {
    double ret;

    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);

    DROPS::FaceCL* face= face_map_[key];

    if (face == NULL) {//non-barycenter
      int ix, iy, iz, it;
      GetNum(p,t,ix,iy,iz,it);
      ret = (*B_Inter_)(ix,0,iz,it);
      //ret= B_Inter_[GetNum(p,t,3)];
    } else {
      int ix1,iy1, iz1,it1;
      GetNum(face->GetVertex(0)->GetCoord(),t,ix1,iy1,iz1,it1);
      int ix2,iy2,iz2,it2;
      GetNum(face->GetVertex(1)->GetCoord(),t,ix2,iy2,iz2,it2);
      int ix3,iy3,iz3,it3;
      GetNum(face->GetVertex(2)->GetCoord(),t,ix3,iy3,iz3,it3);
      iy1 = 0; iy2 = 0; iy3 = 0;
      ret= 1./3.*((*B_Inter_)(ix1,iy1,iz1,it1)+(*B_Inter_)(ix2,iy2,iz2,it2)+(*B_Inter_)(ix3,iy3,iz3,it3));
    }
    return ret;
  };
  //y=0: if dirichlet condition is active
   double GetInterfaceValue( const DROPS::Point3DCL& p, double t) const
  {
    int ix, iy,iz, it;
    GetNum(p,t,ix,iy,iz,it);
    return (*B_Inter_)(ix,0,iz,it);
  };
  //rhs
   double GetSource( const DROPS::Point3DCL& p, double t) const
  {
    double ret;

    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    DROPS::TetraCL* tetra= tetra_map_[key];

    if (tetra == NULL) {//non-barycenter
      int ix1,iy1,iz1,it1;
      GetNum(p,t,ix1,iy1,iz1,it1);
      ret=(*F_)(ix1,iy1,iz1,it1);
    } else {
      int ix1,iy1, iz1,it1;
      GetNum(tetra->GetVertex(0)->GetCoord(),t,ix1,iy1,iz1,it1);
      int ix2,iy2,iz2,it2;
      GetNum(tetra->GetVertex(1)->GetCoord(),t,ix2,iy2,iz2,it2);
      int ix3,iy3,iz3,it3;
      GetNum(tetra->GetVertex(2)->GetCoord(),t,ix3,iy3,iz3,it3);
      int ix4,iy4,iz4,it4;
      GetNum(tetra->GetVertex(3)->GetCoord(),t,ix4,iy4,iz4,it4);
      ret = 0.25*((*F_)(ix1,iy1,iz1,it1)+(*F_)(ix2,iy2,iz2,it2)+
		  (*F_)(ix3,iy3,iz3,it3)+(*F_)(ix4,iy4,iz4,it4));
    }

    return ret;
  };
  //coefficient functions
   double GetDiffusion( const DROPS::Point3DCL& p, double t) const
  {
    double ret;

    d_pair pr= std::make_pair(rnd(p[2]), rnd(p[1]));
    cmp_key key= std::make_pair(rnd(p[0]), pr);
    DROPS::TetraCL* tetra= tetra_map_[key];

    if (tetra == NULL) {//non-barycenter
      int ix1,iy1, iz1,it1;
      GetNum(p,t,ix1,iy1,iz1,it1);
      ret=(*Dw_)(ix1,iy1,iz1,it1)+D_mol_;
    }else {
      int ix1,iy1, iz1,it1;
      GetNum(tetra->GetVertex(0)->GetCoord(),t,ix1,iy1,iz1,it1);
      int ix2,iy2,iz2,it2;
      GetNum(tetra->GetVertex(1)->GetCoord(),t,ix2,iy2,iz2,it2);
      int ix3,iy3,iz3,it3;
      GetNum(tetra->GetVertex(2)->GetCoord(),t,ix3,iy3,iz3,it3);
      int ix4,iy4,iz4,it4;
      GetNum(tetra->GetVertex(3)->GetCoord(),t,ix4,iy4,iz4,it4);
      ret = 0.25*((*Dw_)(ix1,iy1,iz1,it1)+(*Dw_)(ix2,iy2,iz2,it2)+
		  (*Dw_)(ix3,iy3,iz3,it3)+(*Dw_)(ix4,iy4,iz4,it4)) + D_mol_;
    }
    return ret;
  };

  template<class P1EvalT>
    void SetSol3D( const P1EvalT& sol, double t)  //Instationary problem
    {
      bool adjoint = P.get<int>("PoissonCoeff.Adjoint");
      double *out;
      if (adjoint) {
	const int num = (Nt_-rd(t/dt_)-1)*Nxyz_;
	out = C3D_+num;
      } else {
	const int num= (rd(t/dt_)-1)*Nxyz_;  // omit initial time step in output
	out= C3D_+num;               //don't like this way
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

    void SetProductFun(PdeFunction* f1, PdeFunction* f2)
    {
      f1_ = f1;
      f2_ = f2;
    }

  //Check the input matrices
  void Init( const DROPS::ParamCL& P, const PdeFunction* C0, const PdeFunction* B_in, const PdeFunction* F, const PdeFunction* Dw, const PdeFunction* B_Inter, double* c_sol)
  {
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

    // Save the matrix input arguments.
    C0_   = C0;
    B_in_ = B_in;
    F_    = F;
    Dw_    = Dw;
    D_mol_ = P.get<double>("PoissonCoeff.Dmol");
    dt_    = P.get<double>("Time.StepSize");
    B_Inter_= B_Inter;


    // Set the output pointer to the output arguments.
    C3D_ = c_sol;
  }

  void Init(const DROPS::ParamCL& P)
  {
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

    dt_= P.get<double>("Time.StepSize");
  }

};

DROPS::MultiGridCL* PythonConnectCL::MG_= NULL;
PythonConnectCL::FACE_MAP PythonConnectCL::face_map_;
PythonConnectCL::TETRA_MAP PythonConnectCL::tetra_map_;

#endif

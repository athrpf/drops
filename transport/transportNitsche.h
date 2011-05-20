/// \file transportNitsche.h
/// \brief Classes that constitute a 2-phase-transport-problem with Nitsche-XFEM discretization.
/// \author Trung Hieu Nguyen (small fixes: Martin Horsky, Christoph Lehrenfeld), Christoph Lehrenfeld IGPM

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

#ifndef DROPS_TRANSPORTNITSCHE_H
#define DROPS_TRANSPORTNITSCHE_H

#include "levelset/levelset.h"
#include "levelset/mgobserve.h"
#include "misc/problem.h"
#include "num/fe.h"
#include "num/discretize.h"
#include "num/solver.h"
#include "num/bndData.h"
#include "stokes/instatstokes2phase.h"
#include "misc/xfem.h"
#include <iostream>
#include <numeric>
#include <cstring>

namespace DROPS
{
  
class VelocityContainer
{
  typedef BndDataCL<Point3DCL>                              VelBndDataT;
  typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;  
  private:
    VecDescCL *v_;   
    const VelBndDataT*    Bnd_v_;     
    MultiGridCL*   MG_;  
    instat_vector_fun_ptr vfptr_; 
    const_DiscVelSolCL * asp2;
  public:
    VelocityContainer(VecDescCL & v,const VelBndDataT& Bnd_v,MultiGridCL& MG):v_(&v),Bnd_v_(&Bnd_v), MG_(&MG),vfptr_(0)
    {
      asp2 = new const_DiscVelSolCL( v_, Bnd_v_, MG_);
    };
    
    VelocityContainer(instat_vector_fun_ptr v):v_(0),Bnd_v_(0),MG_(0),vfptr_(v),asp2(0){};
    
    ~VelocityContainer()
    {
      if (asp2) delete asp2;
    }
    
    const_DiscVelSolCL & GetVelocityAsP2() const
        { 
          if (!(v_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return *asp2; 
        }
        
    const_DiscVelSolCL GetVelocityAsP2(const VecDescCL& vel) const
        { 
          if (!(v_ && Bnd_v_))
            throw DROPSErrCL("velocity not prescribed as a const_DiscVelSolCL");
          return const_DiscVelSolCL( &vel, Bnd_v_, MG_); 
        }

    instat_vector_fun_ptr GetVelocityAsFunctionPointer() const
    {
      if (!vfptr_)
        throw DROPSErrCL("velocity not prescribed as a function(pointer)");
      return vfptr_;
    }
    
    bool VelocityAsP2() const {
      return (v_ && Bnd_v_);
    }
     
    bool VelocityAsFunctionPointer() const {
      return (vfptr_);
    }
 
};
  
class TransportP1XCL
{
  public:
    typedef BndDataCL<>                                       BndDataT;
    typedef BndDataCL<Point3DCL>                              VelBndDataT;
    typedef P1EvalCL<double, const BndDataT, VecDescCL>       DiscSolCL;
    typedef P1EvalCL<double, const BndDataT, const VecDescCL> const_DiscSolCL;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;

    double oldt_, t_;
    MLIdxDescCL idx, oldidx;
    VecDescCL c,           ///< concentration
              ct, oldct;   ///< transformed concentration
    MLMatDescCL A,         ///< diffusion matrix
                M,         ///< mass matrix
                C,        ///< convection matrix
                NA;       ///< Nitsche term 
    VecDescCL cplA, cplM, cplC, b; 
    
  private:
    MLMatrixCL     L_;              ///< sum of matrices
    MultiGridCL&   MG_;
    BndDataT&      Bnd_, Bndt_;     ///< Boundary condition for the concentration
    VelocityContainer v_;
    LsetBndDataCL& Bnd_ls_;          ///< Boundary condition for the level set function
    double         theta_,          ///< time scheme parameter
                   dt_,             ///< time step 
                   lambda_,         ///< Nitsche's parameter
                   D_[2],           ///< diffusion constants
                   H_; 
    VecDescCL &lset_,  &oldlset_;  ///< levelset at current time step and previous step
    GSPcCL                  pc_;
    GMResSolverCL<GSPcCL>   gm_;
    instat_scalar_fun_ptr f_;         ///<source term
    instat_scalar_fun_ptr c_;        ///<mass/reaction term
    double omit_bound_;              ///discard criteria for XFEM-Basis
    double sdstab_;
    void SetupInstatSystem(MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, VecDescCL*,
        IdxDescCL&, const double) const;
    void SetupNitscheSystem( MatrixCL&, IdxDescCL&)const;
    void SetupNitscheSystem( MatrixCL&, IdxDescCL&, bool)const;
    void SetupInstatMixedMassMatrix( MatrixCL&, VecDescCL*, IdxDescCL&, IdxDescCL&, const double) const;
    void SetupInstatMixedSystem(MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, MatrixCL&, VecDescCL*, VecDescCL*,
        IdxDescCL&, IdxDescCL&, const double) const;
    void SetupMixedNitscheSystem( MatrixCL&, IdxDescCL&, IdxDescCL&, const double) const;

  public:
/*    TransportP1XCL( MultiGridCL& mg, BndDataT& Bnd, BndDataT& Bndt, const VelBndDataT& Bnd_v, LsetBndDataCL& Bnd_ls,
        double theta, double D[2], double H, VecDescCL* v, VecDescCL* oldv, VecDescCL& lset, VecDescCL& oldlset,
        double dt, instat_scalar_fun_ptr rhs=0, int iter= 1000, double tol= 1e-7,
        double XFEMstab=0.1, double lambda=0.)
        : oldt_(v->t), t_( v->t), 
        idx( P1X_FE, 1, Bndt, 0, XFEMstab), oldidx( P1X_FE, 1, Bndt, 0, XFEMstab), 
        MG_( mg), Bnd_( Bnd), Bndt_( Bndt), Bnd_v_( &Bnd_v), Bnd_ls_(Bnd_ls), 
        theta_( theta), dt_( dt), 
        lambda_(lambda), H_( H), v_( v), oldv_(oldv), 
        lset_( lset), oldlset_(oldlset), 
        gm_( pc_, 20, iter, tol, false, false, RightPreconditioning),
        f_(rhs) 
    {
        std::memcpy( D_, D, 2*sizeof( double));
        if (theta_!=1.0) std::cerr<< "TransportP1XCL::TransportP1XCL: Sorry. Cannot deal with theta != 1.0. Overwriting: Using implicit Euler now! \n";
    }
*/
    TransportP1XCL( MultiGridCL& mg, BndDataT& Bnd, BndDataT& Bndt, VelocityContainer& v, LsetBndDataCL& Bnd_ls,
        VecDescCL& lset, VecDescCL& oldlset,
        DROPS::ParamMesszelleNsCL& C, double initialtime=0, instat_scalar_fun_ptr reac=0, instat_scalar_fun_ptr rhs=0)
        : oldt_(initialtime), t_( initialtime), 
        idx( P1X_FE, 1, Bndt, mg.GetBnd().GetMatchFun(), C.trp_NitscheXFEMStab), oldidx( P1X_FE, 1, Bndt, mg.GetBnd().GetMatchFun(), C.trp_NitscheXFEMStab), 
        MG_( mg), Bnd_( Bnd), Bndt_( Bndt), v_ (v), Bnd_ls_(Bnd_ls), 
        theta_( C.trp_Theta), dt_( C.tm_StepSize), 
        lambda_(C.trp_NitschePenalty), H_( C.trp_HNeg/C.trp_HPos),
        lset_( lset), oldlset_(oldlset), 
        gm_( pc_, 20, C.trp_Iter, C.trp_Tol, false, false, RightPreconditioning),
        f_(rhs), c_(reac), omit_bound_( C.trp_NitscheXFEMStab), sdstab_(C.trp_SDStabilization)
    {
        double D[2] = {C.trp_DiffPos, C.trp_DiffNeg};
        std::memcpy( D_, D, 2*sizeof( double));
        if (theta_!=1.0) std::cerr<< "TransportP1XCL::TransportP1XCL: Sorry. Cannot deal with theta != 1.0. Overwriting: Using implicit Euler now! \n";
    }

    const VelocityContainer & GetVelocity() const {return v_;}

    MultiGridCL& GetMG()               { return MG_; }
    const MultiGridCL& GetMG() const   { return MG_; }
    GMResSolverCL<GSPcCL>& GetSolver() { return gm_; }
    const BndDataT& GetBndData() const { return Bndt_; }
    
    /// \brief Only used for XFEM
    void UpdateXNumbering( const VecDescCL& lset, bool= false, bool NewSolution=true) {
         if (NewSolution)
            idx.UpdateXNumbering( MG_, lset, Bnd_ls_);
         else
            oldidx.UpdateXNumbering( MG_, lset, Bnd_ls_);
    }

    void CreateNumbering( Uint level, MLIdxDescCL* idx1, MLIdxDescCL *oldidx1, VecDescCL& lset1, VecDescCL& oldlset1) {
         idx1->CreateNumbering( level, MG_, &lset1, &Bnd_ls_);
         oldidx1->CreateNumbering( level, MG_, &oldlset1, &Bnd_ls_);
    }
    void DeleteNumbering( MLIdxDescCL* idx1) { idx1->DeleteNumbering( MG_); }
    /// initialize transformed concentration function
    void Init( instat_scalar_fun_ptr, instat_scalar_fun_ptr, double t=0.0);
    /// Set one P1X Function as the other with according scaling parameters in the positive and negative domain parts
    void TransformWithScaling(const VecDescCL& concin, VecDescCL& concout, double scalingp, double scalingn);

    void SetTimeStep( double dt, double theta=-1);
    void SetTwoStepIdx();
    void SetNewIdx();
    void SetupInstatSystem( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, VecDescCL&, const double) const;
    void SetupNitscheSystem( MLMatDescCL&) const;
    void SetupNitscheSystem( MLMatDescCL&, bool) const;
    void SetupInstatMixedMassMatrix( MLMatDescCL&, VecDescCL&, const double) const;
    void SetupInstatMixedSystem( MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, MLMatDescCL&, VecDescCL&, VecDescCL&, const double) const;
    void SetupMixedNitscheSystem( MLMatDescCL&, const double) const;

    /// perform one time step
    void DoStep( double new_t);
    void SetRHS( instat_scalar_fun_ptr rhs) {f_= rhs;}
    const_DiscSolCL GetSolution() const
        { return const_DiscSolCL( &c, &Bnd_, &MG_); }
    const_DiscSolCL GetSolution( const VecDescCL& Myc, bool Is_ct) const
        { return (Is_ct) ? const_DiscSolCL( &Myc, &Bndt_, &MG_)
                         : const_DiscSolCL( &Myc, &Bnd_, &MG_); }
                         
 
   
    /// \name For internal use only
    /// The following member functions are added to enable an easier implementation
    /// of the locling navstokes-levelset. They should not be called by a common user.
    /// Use LevelsetP2CL::DoStep() instead.
    ///@{
    void InitStep (VectorCL&);
    void DoStep (const VectorCL&);
    void CommitStep ();
    ///@}
    /// Get FE type for pressure space
    VecDescCL& GetLevelset ()  {return lset_;}
    const VecDescCL& GetLevelset ()  const {return lset_;}
    LsetBndDataCL& GetLevelsetBnd ()  {return Bnd_ls_;}
    const LsetBndDataCL& GetLevelsetBnd ()  const {return Bnd_ls_;}
    VecDescCL& GetOldLevelset ()  {return oldlset_;}
    const VecDescCL& GetOldLevelset ()  const {return oldlset_;}
    //@}
    void GetSolutionOnPart( VecDescCL& ct_part, bool posPart, bool Is_ct );
    double Interface_L2error() const;
    void SetupInstatRhs( VecDescCL & b,  const double time) const;
    double MeanDropConcentration();
    double GetHenry(bool positive = false){
      if (positive)
        return 1.0;
      else
        return H_;
    }
    double GetXFEMOmitBound(){ return omit_bound_; }
    double CheckSolution(instat_scalar_fun_ptr Lsgn, instat_scalar_fun_ptr Lsgp, double time);
};

/// \brief Observes the MultiGridCL-changes by AdapTriangCL to repair the function c.ct.
///
/// The actual work is done in post_refine().

class TransportXRepairCL : public MGObserverCL
{
  private:
    TransportP1XCL& c_;
    std::auto_ptr<P1XRepairCL> oldp1xrepair_;
    Uint mylvl;
  public:
    TransportXRepairCL (TransportP1XCL& c, Uint amylvl)
        : c_( c), mylvl(amylvl) {}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  (); 
    void post_refine_sequence (); 
    const IdxDescCL* GetIdxDesc() const;
};

class VelTranspRepairCL : public MGObserverCL
{
  public:
    typedef BndDataCL<Point3DCL>                                          VelBndDataT;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL>    const_DiscVelSolCL;
  private:
    MultiGridCL& mg_;
    VecDescCL& v_;
    VecDescCL& lset_;
    const VelBndDataT& Bnd_v_;
    LsetBndDataCL& Bnd_ls_;
    IdxDescCL& vidx_;
    double& time_;

  public:
    VelTranspRepairCL (VecDescCL& v, MultiGridCL& mg, const VelBndDataT& Bnd_v, IdxDescCL& vidx, VecDescCL& lset, LsetBndDataCL& Bnd_ls, double t )
        :  mg_(mg), v_(v), lset_(lset), Bnd_v_(Bnd_v), Bnd_ls_(Bnd_ls), vidx_(vidx) , time_(t){}

    void pre_refine  () {}
    void post_refine ();

    void pre_refine_sequence  () {}
    void post_refine_sequence () {}
    void SetTime(double t) {time_=t;}
    const IdxDescCL* GetIdxDesc() const;
};


/// \todo: think about where to put this in a more systematical way

class P1FEGridfunctions{
  private:
    LocalP1CL<> p1[4];
    LocalP2CL<> p2[4];
    Quad3CL<> q3_p[4];
    Quad5CL<> q5_p[4];
    LocalP2CL<> pipj[4][4];
  public:
    static void SetupPiPj(LocalP2CL<>pipj[4][4])
    {    
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < i; ++j) {
                pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
                pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
            }
            pipj[i][i][i]= 1.;
            for (int vert= 0; vert < 3; ++vert)
                pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
        }
    }

    P1FEGridfunctions(){
      for(int i= 0; i < 4; ++i) {
          p1[i][i]=1.;
          q3_p[i].assign(p1[i]);
          q5_p[i].assign(p1[i]);
          p2[i].assign(p1[i]);
       }
      SetupPiPj(pipj);  
    }
    
    LocalP1CL<>& GetShapeAsLocalP1(int i) {return p1[i];}
    LocalP2CL<>& GetShapeAsLocalP2(int i) {return p2[i];}
    LocalP2CL<>& GetProductShapeAsLocalP2(int i, int j) {return pipj[i][j];}
    Quad3CL<>& GetShapeAsQuad3(int i) {return q3_p[i];}
    Quad5CL<>& GetShapeAsQuad5(int i) {return q5_p[i];}
};

class LocalConvDiffReacCoefficients; //forward declaration

class TransformedP1FiniteElement{
  private:
    SMatrixCL<3,4> G;
    SMatrixCL<4,4> GTG;
    double det, absdet;
    double vol;
  protected:
    TetraCL* tet;
    bool has_trafo_base;
    bool has_Gram;
    bool oninterface;
    P1FEGridfunctions& p1fegfs;
    BaryCoordCL* nodes;
    Quad3CL<> q3_baseshape[4];
  public:
    TransformedP1FiniteElement(P1FEGridfunctions& ap1fegfs, TetraCL* atet = NULL):tet(atet), p1fegfs(ap1fegfs){
      has_trafo_base=false;    
      has_Gram=false;    
      oninterface = false; 
      nodes = NULL;
    }
    
    ~TransformedP1FiniteElement(){
      if (nodes) delete nodes;
    }
    
    P1FEGridfunctions& GetGridfunctions(){
      return p1fegfs;
    }
    
    void SetTetra(TetraCL& atet){
      tet = &atet;
      has_trafo_base=false;    
      has_Gram=false;    
      oninterface = false;    
    }
    
    virtual void SetLocal(TetraCL& atet, LocalConvDiffReacCoefficients& , bool ){    
      SetTetra(atet);
    }

    void SetSubTetra(const SArrayCL<BaryCoordCL,4>& cutT){
      vol = VolFrac(cutT) * absdet * 1.0/6.0;
      if (nodes) delete nodes;
      nodes = Quad3CL<>::TransformNodes(cutT);
      oninterface = true;
      
      for (int i = 0; i < 4; i ++){
        q3_baseshape[i] = Quad3CL<>(p1fegfs.GetShapeAsLocalP2(i), nodes);      
      }      
    }

    TetraCL& GetTetra() const{
      if(tet == NULL) throw DROPSErrCL("TransformedP1FiniteElement::GetTetra - No TetraCL object given!");
      return *tet;
    }
    
    void CalcTrafoBase(){
      P1DiscCL::GetGradients(G, det, GetTetra());
      absdet = std::fabs( det);
      vol = absdet * 1.0/6.0;
      has_trafo_base = true;
    }
    
    BaryCoordCL* GetNodes() const{
      if (!oninterface)
        throw DROPSErrCL("GetNodes should only be called if a tetra is subdivided!");
      else
        return nodes;
    }
    
    double GetDeterminant(){
      if (!has_trafo_base) CalcTrafoBase();
      return det;
    }

    double GetAbsDeterminant(){
      if (!has_trafo_base) CalcTrafoBase();
      return absdet;
    }

    double GetVolume(){
      if (!has_trafo_base) CalcTrafoBase();
      return vol;
    }

    SMatrixCL<3,4> & GetDShape(){
      if (!has_trafo_base) CalcTrafoBase();
      return G;
    }

    SMatrixCL<4,4> & GetGramShape(){
      if (!has_Gram){
        if (!has_trafo_base) CalcTrafoBase();
        GTG = GramMatrix(G);
        has_Gram = true;
      } 
      return GTG;
    }

    Quad3CL<> & GetBaseShapeAsQuad3CL(int i) {
      if (!oninterface)
        return p1fegfs.GetShapeAsQuad3(i);
      else
        return q3_baseshape[i];
    }

    virtual bool stabilized() const{
      return false;
    }
  
};

class GlobalConvDiffReacCoefficients{
  friend class LocalConvDiffReacCoefficients;
  private:
    typedef BndDataCL<Point3DCL> VelBndDataT;
    typedef P2EvalCL<SVectorCL<3>, const VelBndDataT, const VecDescCL> const_DiscVelSolCL;  
    double D_[2];
    double H_;
    double time;    
    const VelocityContainer& vel;
    instat_scalar_fun_ptr source;
    instat_scalar_fun_ptr mass;
  public:

    GlobalConvDiffReacCoefficients(const double D[2], double H, const VelocityContainer& u, instat_scalar_fun_ptr c, instat_scalar_fun_ptr f, double atime)
      : H_(H),time(atime),vel(u),source(f),mass(c){
      D_[0]=D[0];
      D_[1]=D[1];
    }

    double GetHenryWeighting(bool pospart){
      if (pospart)
        return 1.0;
      else
        return H_;
    }

    double GetDiffusionCoef(bool pospart){
      if (pospart)
        return D_[0];
      else
        return D_[1];
    }
	
	instat_scalar_fun_ptr GetSourceAsFuntionPointer(){return source;}
	instat_scalar_fun_ptr GetMassAsFuntionPointer(){return mass;}
};

class LocalConvDiffReacCoefficients{
  private:
    GlobalConvDiffReacCoefficients& gcdcoefs;  
    Quad5CL<> q5_source;
    Quad3CL<> q3_source;
    Quad5CL<> q5_mass;
    Quad3CL<> q3_mass;
    Quad3CL<Point3DCL> *q3_velocity;
    LocalP2CL<Point3DCL> *lp2_velocity;  
  public:
    LocalConvDiffReacCoefficients(GlobalConvDiffReacCoefficients& agcdcoefs, TetraCL& tet):gcdcoefs(agcdcoefs), 
        q5_source(tet,gcdcoefs.source,gcdcoefs.time), q3_source(tet,gcdcoefs.source,gcdcoefs.time),
        q5_mass(tet,gcdcoefs.mass,gcdcoefs.time), q3_mass(tet,gcdcoefs.mass,gcdcoefs.time)
    {
      if (gcdcoefs.vel.VelocityAsP2()){
        q3_velocity = new Quad3CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsP2());
        lp2_velocity = new LocalP2CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsP2());
      }
      else
      {
        q3_velocity = new Quad3CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsFunctionPointer(),gcdcoefs.time);
        lp2_velocity = new LocalP2CL<Point3DCL>(tet, gcdcoefs.vel.GetVelocityAsFunctionPointer(),gcdcoefs.time);
      }
    }
    
    ~LocalConvDiffReacCoefficients(){
      delete q3_velocity;
      delete lp2_velocity;
    }
    double GetHenryWeighting(bool pospart){
      return gcdcoefs.GetHenryWeighting(pospart);
    }

    double GetDiffusionCoef(bool pospart){
      return gcdcoefs.GetDiffusionCoef(pospart);
    }    
	
	instat_scalar_fun_ptr GetSourceAsFuntionPointer(){
	  return gcdcoefs.GetSourceAsFuntionPointer();
	}
	double GetTime(){return gcdcoefs.time;}
	Quad5CL<>& GetSourceAsQuad5(){return q5_source;}
	Quad3CL<>& GetSourceAsQuad3(){return q3_source;}
	Quad5CL<>& GetReactionAsQuad5(){return q5_mass;}
	Quad3CL<>& GetReactionAsQuad3(){return q3_mass;}
	Quad3CL<Point3DCL>& GetVelocityAsQuad3(){return *q3_velocity;}
	LocalP2CL<Point3DCL>& GetVelocityAsLocalP2(){return *lp2_velocity;}

};


class StabilizedTransformedP1FiniteElement : public TransformedP1FiniteElement{
  protected:
    LocalConvDiffReacCoefficients* lcdcoefs;
    double h;
    
    double delta_T; ///element-wise constant stabilization parameter
    
    Quad3CL<> q3_stabbgradv[4];
    Quad3CL<> q3_testshape[4];
//    Quad3CL<> q3_baseshape[4];
//    Quad5CL<> q5_stabbgradv[4];
    double stabfactor_;
  public:
    typedef TransformedP1FiniteElement base_;  
    StabilizedTransformedP1FiniteElement(P1FEGridfunctions& ap1fegfs, double stabfactor = 1.0, LocalConvDiffReacCoefficients* alcdcoefs= NULL, TetraCL* atet = NULL): base_(ap1fegfs,atet),lcdcoefs(alcdcoefs){
      stabfactor_ = stabfactor;
      std::cout << "stabfactor_ = " << stabfactor_ << std::endl;
      has_trafo_base=false;    
      has_Gram=false;  
      delta_T = 0.0;
    }
    
    double GetStabilizationParameter() const{
      return delta_T;
    }
    
    double stabfunc(double diffusion, double velocity, double meshsize) const {
      double meshPeclet = meshsize * 0.5 * velocity / diffusion;
      double res = 0;
      if (meshPeclet > 1.0){
        res = (1 - 1.0/meshPeclet) * meshsize/(2*velocity);
        //std::cout << "meshPeclet = " << meshPeclet << std::endl;
      }
      //return 0.;
      return res;
    }
    
    void CalcStabilizationOnePhase(bool pospart){
      CalcTrafoBase(); //s.t. gradients, etc.. are known.
      //l2-average of velocity
      double ul2 = 1.0/lcdcoefs->GetHenryWeighting(pospart)*std::sqrt( (Quad3CL<>(dot( lcdcoefs->GetVelocityAsQuad3(), lcdcoefs->GetVelocityAsQuad3()))).quad(6.0));
      double h = cbrt(6.0*GetVolume()); 
      delta_T = stabfactor_ * stabfunc(lcdcoefs->GetDiffusionCoef(pospart),ul2,h);
      //ATTENTION: convdiff-debug
#ifdef ONLY_OUTERCONVECTION      
      if (!pospart) delta_T = 0;
#endif
      Point3DCL stabtimesgradv;
      SMatrixCL<3,4> & dshape = GetDShape();
      for (int i = 0; i < 4; i ++){
        for (int j = 0; j < 3; j ++)
          stabtimesgradv[j] = dshape(j,i);
        q3_stabbgradv[i] = dot(stabtimesgradv, lcdcoefs->GetVelocityAsQuad3());
        q3_testshape[i] = delta_T * q3_stabbgradv[i] + base_::GetGridfunctions().GetShapeAsQuad3(i);
      }
    }
    
    void CalcStabilizationTwoPhase(bool pospart){
      //l2-average of velocity
      Quad3CL<Point3DCL> q3_n_vel(lcdcoefs->GetVelocityAsLocalP2(), nodes);      
      double ul2 = 1.0/lcdcoefs->GetHenryWeighting(pospart)*std::sqrt( (Quad3CL<>(dot( q3_n_vel, q3_n_vel))).quad(6.0));
      double h = cbrt(6.0*GetVolume());
      delta_T = stabfactor_ * stabfunc(lcdcoefs->GetDiffusionCoef(pospart),ul2,h);

      //ATTENTION: convdiff-debug
      //if (!pospart) delta_T = 0;
      Point3DCL stabtimesgradv;
      SMatrixCL<3,4> & dshape = GetDShape();
      for (int i = 0; i < 4; i ++){
        Quad3CL<> q3_n_test(base_::GetGridfunctions().GetShapeAsLocalP2(i), nodes);      
        for (int j = 0; j < 3; j ++)
          stabtimesgradv[j] = dshape(j,i);
        q3_stabbgradv[i] = dot(stabtimesgradv, q3_n_vel);
        q3_baseshape[i] = q3_n_test;
        q3_testshape[i] = delta_T * q3_stabbgradv[i] + q3_baseshape[i];
      }
    }
     
    P1FEGridfunctions& GetGridfunctions() {
      return base_::GetGridfunctions();
    }
    
    virtual void SetLocal(TetraCL& atet, LocalConvDiffReacCoefficients& alcdcoefs, bool pospart){
      base_::SetTetra(atet);
      lcdcoefs = &alcdcoefs;
      oninterface = false;
      CalcStabilizationOnePhase(pospart);
    }
    
    void CalcStabilization(bool pospart){
      if(oninterface)
        CalcStabilizationTwoPhase(pospart);
      else 
        CalcStabilizationOnePhase(pospart);
    };

    Quad3CL<> & GetStabTestShapeAsQuad3CL(int i) {
      return q3_stabbgradv[i];
    }

    Quad3CL<> & GetTestShapeAsQuad3CL(int i) {
      return q3_testshape[i];
    }
    

    virtual bool stabilized() const{
//      return false;
//      return true;
//
      return (delta_T > 0);
    }

};





typedef double Elmat4x4[4][4];
typedef double Elvec4[4];


class ConvDiffElementVectors{
  public:

    Elvec4 f;
    Elvec4 f_p;
    Elvec4 f_n;

    void ResetUnsigned(){
      std::memset( f,0, 4*sizeof(double));
    }
    
    void ResetSigned(){
      std::memset( f_p,0, 4*sizeof(double));
      std::memset( f_n,0, 4*sizeof(double));  
    }

    void ResetAll(){
      ResetUnsigned();
      ResetSigned();
    }
    
    ConvDiffElementVectors(){
      ResetAll();
    }
};


class ConvDiffElementMatrices{
  public:
    Elmat4x4 A;
    Elmat4x4 M;
    Elmat4x4 Mr;
    Elmat4x4 C;
    Elmat4x4 A_p;
    Elmat4x4 M_p;
    Elmat4x4 Mr_p;
    Elmat4x4 C_p;
    Elmat4x4 A_n;
    Elmat4x4 M_n;
    Elmat4x4 Mr_n;
    Elmat4x4 C_n;

    void ResetUnsigned(){
      std::memset( A,0, 4*4*sizeof(double));
      std::memset( M,0, 4*4*sizeof(double));
      std::memset( Mr,0, 4*4*sizeof(double));
      std::memset( C,0, 4*4*sizeof(double));
    }
    
    void ResetSigned(){
      std::memset( A_p,0, 4*4*sizeof(double));
      std::memset( M_p,0, 4*4*sizeof(double));
      std::memset( Mr_p,0, 4*4*sizeof(double));
      std::memset( C_p,0, 4*4*sizeof(double));
      std::memset( A_n,0, 4*4*sizeof(double));
      std::memset( M_n,0, 4*4*sizeof(double));
      std::memset( Mr_n,0, 4*4*sizeof(double));
      std::memset( C_n,0, 4*4*sizeof(double));
    }

    void ResetAll(){
      ResetUnsigned();
      ResetSigned();
    }
    
    void SetUnsignedAsSumOfSigned(){
      for(int i= 0; i < 4; ++i){
        for(int j= 0; j < 4; ++j){
          M[i][j]= M_n[i][j] + M_p[i][j];
          A[i][j]= A_n[i][j] + A_p[i][j];
          C[i][j]= C_n[i][j] + C_p[i][j];
          Mr[i][j]= Mr_n[i][j] + Mr_p[i][j];
        }
      }
    }
    
    ConvDiffElementMatrices(){
      ResetAll();
    }
};



} // end of namespace DROPS

#endif

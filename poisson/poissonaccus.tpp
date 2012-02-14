/// \file poissonaccus.tpp
/// \brief classes that constitute the integrals of the finite element formulation for the poisson-problem
/// \author LNM RWTH Aachen: Christoph Lehrenfeld, Liang Zhang

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
namespace DROPS
{

/// - Base Accumulator for scalar P1 finite elements

template<class Coeff,template <class T=double> class QuadCL>
class Accumulator_P1CL : public TetraAccumulatorCL
{
    protected:   
    const MultiGridCL& MG_;
    const BndDataCL<> * BndData_; 
    MatrixCL* Amat_; 
    VecDescCL* b_; 
    IdxDescCL& RowIdx_; 
    IdxDescCL& ColIdx_; 
    
    MatrixBuilderCL * A_;    
    
    //local informations
    
    // - sharable (for future changes)
    Point3DCL G[4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];
    LocalP1CL<double> phi[4];
    QuadCL<> phiQuad[4];
    // - not sharable (for future changes)
    double coup[4][4];
    QuadCL<> U_Grad[4];
    
    const Uint lvl;
    const Uint idx;
    
    const double t;
    void update_global_matrix();
    void update_coupling(const TetraCL& sit);
    
    public: 
        Accumulator_P1CL (const MultiGridCL& MG, const BndDataCL<> * BndData, MatrixCL* Amat, VecDescCL* b, 
                   IdxDescCL& RowIdx, IdxDescCL& ColIdx, const double t_=0);

    ///\brief Initializes matrix-builders and load-vectors
    void begin_accumulation ();
    ///\brief Builds the matrices
    void finalize_accumulation();

    virtual void visit (const TetraCL&) 
        { throw DROPSErrCL("BaseClass Accumulator_P1CL::visit called - this should not happen!");};

    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new Accumulator_P1CL ( *this); }
        
    
};

template<class Coeff,template <class T=double> class QuadCL>
void Accumulator_P1CL<Coeff,QuadCL>::update_global_matrix()
{
    for(int i=0; i<4; ++i)          // assemble row i
        if (UnknownIdx[i]!= NoIdx)  // vertex i is not on a Dirichlet boundary
        {
            for(int j=0; j<4;++j)
            {
                if (UnknownIdx[j]!= NoIdx) // vertex j is not on a Dirichlet boundary
                {
                    (*A_)(UnknownIdx[i], UnknownIdx[j])+=coup[i][j];    //DiffusiconCoeff*A
                }
            }
        } 
}

template<class Coeff,template <class T=double> class QuadCL>
void Accumulator_P1CL<Coeff,QuadCL>::update_coupling(const TetraCL& sit)
{
    for(int i=0; i<4; ++i)          // assemble row i
        if (UnknownIdx[i]!= NoIdx)  // vertex i is not on a Dirichlet boundary
        {
            for(int j=0; j<4;++j)
            {
                if (UnknownIdx[j]== NoIdx) // vertex j is on a Dirichlet boundary
                {
                    b_->Data[UnknownIdx[i]]-= coup[i][j] * BndData_->GetDirBndValue(*sit.GetVertex(j), t);
                }
            }
        } 
}

template<class Coeff,template <class T=double> class QuadCL>
Accumulator_P1CL<Coeff,QuadCL>::Accumulator_P1CL(const MultiGridCL& MG, const BndDataCL<> * BndData, MatrixCL* Amat, VecDescCL* b, 
                   IdxDescCL& RowIdx, IdxDescCL& ColIdx, const double t_):                  
                   MG_(MG), BndData_(BndData), Amat_(Amat), b_(b), RowIdx_(RowIdx), ColIdx_(ColIdx), 
                   A_(0),
                   lvl(RowIdx.TriangLevel()),
                   idx(RowIdx.GetIdx()), t(t_)
{
    for(int i=0; i<4; i++)
    {
      phi[i][i]=1.;
      phiQuad[i].assign(phi[i]);
    }
   
}
template<class Coeff,template <class T=double> class QuadCL>
void Accumulator_P1CL<Coeff,QuadCL>::begin_accumulation ()
{
    if (b_ != 0) b_->Clear( t);
    if (Amat_)    
        A_ = new MatrixBuilderCL( Amat_, RowIdx_.NumUnknowns(), ColIdx_.NumUnknowns());
        
}
template<class Coeff,template <class T=double> class QuadCL>
void Accumulator_P1CL<Coeff,QuadCL>::finalize_accumulation ()
{
    if (A_ != 0){
        A_->Build();
        delete A_;   
    }     
    
}

/// - Accumulators for source, stiffness, mass and convection:


// source matrix: \int_{\Omega} f v \, dx
template<class Coeff,template <class T=double> class QuadCL>
class SourceAccumulator_P1CL : public Accumulator_P1CL<Coeff,QuadCL>
{
    protected:
    typedef Accumulator_P1CL<Coeff,QuadCL> base_;
    using                           base_::MG_;
    using                           base_::BndData_; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::A_;              //Stiffnesss matrix
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::phiQuad;
    using                           base_::coup;
    using                           base_::U_Grad;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::t;
    SUPGCL& supg_;
    bool    ALE_;
    QuadCL<> rhs;
    public:
    SourceAccumulator_P1CL(const MultiGridCL& MG, const BndDataCL<> * BndData, VecDescCL* b, 
                   IdxDescCL& RowIdx, SUPGCL& supg, bool ALE, const double t_)
                   :Accumulator_P1CL<Coeff,QuadCL>(MG,BndData,0,b,RowIdx,RowIdx,t_),supg_(supg), ALE_(ALE){}
    void local_setup (const TetraCL& sit);
    void update_rhsintegrals(const TetraCL& sit);
    void visit (const TetraCL& sit);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new SourceAccumulator_P1CL ( *this); }
    
};



template<class Coeff,template <class T=double> class QuadCL>
void SourceAccumulator_P1CL<Coeff,QuadCL>::visit (const TetraCL& sit)
{
  if (b_ != 0 && BndData_ != 0){      
    local_setup(sit);
    update_rhsintegrals(sit);
  }
}

template<class Coeff,template <class T=double> class QuadCL>
void SourceAccumulator_P1CL<Coeff,QuadCL>::local_setup (const TetraCL& sit)
{
    rhs.assign( sit, Coeff::f, t);
    
    for(int i=0; i<4; ++i)
    {
      UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx) : NoIdx;      
    }    
    
    P1DiscCL::GetGradients(G,det,sit);
    absdet= std::fabs(det);    
    
    if(supg_.GetSUPG())
    {    
        instat_vector_fun_ptr vel;
        if(ALE_)
        vel= Coeff::ALEVelocity;
        else
        vel= Coeff::Vel;
        QuadCL<Point3DCL> u(sit,vel,t);
        for(int i=0; i<4; ++i)
            U_Grad[i]=dot( u, QuadCL<Point3DCL>( G[i]));
    }    
}

template<class Coeff,template <class T=double> class QuadCL>
void SourceAccumulator_P1CL<Coeff,QuadCL>::update_rhsintegrals(const TetraCL& sit)
{
    for(int i=0; i<4; ++i)    // assemble row i
        if (UnknownIdx[i]!= NoIdx)  // vertex i is not on a Dirichlet boundary
        {
            QuadCL<double> fp1(rhs*phiQuad[i]);
            b_->Data[UnknownIdx[i]]+= fp1.quad(absdet);
            if (supg_.GetSUPG()) {
                instat_vector_fun_ptr vel;
                if(ALE_)
                vel= Coeff::ALEVelocity;
                else
                vel= Coeff::Vel;
                QuadCL<double> f_SD( rhs*U_Grad[i] );    //SUPG for source term
                b_->Data[UnknownIdx[i]]+= f_SD.quad(absdet)*supg_.Sta_Coeff( vel(GetBaryCenter(sit), t), Coeff::alpha(GetBaryCenter(sit), t));
            }
            if ( BndData_!=0 && BndData_->IsOnNatBnd(*sit.GetVertex(i)) )
                for (int f=0; f < 3; ++f)
                    if ( sit.IsBndSeg(FaceOfVert(i, f)) )
                        b_->Data[UnknownIdx[i]]+= P1DiscCL::Quad2D(sit, FaceOfVert(i, f), BndData_->GetBndSeg(sit.GetBndIdx(FaceOfVert(i,f))).GetBndFun(), i );
        } 
    
}
    


// stiffness matrix: \int_{\Omega} \alpha \nabla u \nabla v \, dx
template<class Coeff,template <class T=double> class QuadCL>
class StiffnessAccumulator_P1CL : public Accumulator_P1CL<Coeff,QuadCL>
{
    protected:
    typedef Accumulator_P1CL<Coeff,QuadCL> base_;
    using                           base_::MG_;
    using                           base_::BndData_; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::A_;              //Stiffnesss matrix
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::phiQuad;
    using                           base_::coup;
    using                           base_::U_Grad;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::t;
    SUPGCL& supg_;
    bool   ALE_;
    public:
    StiffnessAccumulator_P1CL(const MultiGridCL& MG, const BndDataCL<> * BndData, MatrixCL* Amat, VecDescCL* b, 
                   IdxDescCL& RowIdx, IdxDescCL& ColIdx, SUPGCL& supg, bool ALE, const double t_)
                   :Accumulator_P1CL<Coeff,QuadCL>(MG,BndData,Amat,b,RowIdx,ColIdx,t_),supg_(supg), ALE_(ALE){}
    void local_setup (const TetraCL& sit);
    void visit (const TetraCL& sit);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new StiffnessAccumulator_P1CL ( *this); }
    
};


template<class Coeff,template <class T=double> class QuadCL>
void StiffnessAccumulator_P1CL<Coeff,QuadCL>::visit (const TetraCL& sit)
{
  local_setup(sit);
  if (A_ != 0)
    base_::update_global_matrix();
  if (b_ != 0 && BndData_ != 0)      
    base_::update_coupling(sit);
}

template<class Coeff,template <class T=double> class QuadCL>
void StiffnessAccumulator_P1CL<Coeff,QuadCL>::local_setup (const TetraCL& sit)
{    
    Quad2CL<> quad_a;
    P1DiscCL::GetGradients(G,det,sit);
    absdet= std::fabs(det);
    quad_a.assign( sit, Coeff::alpha, 0.0);                  //for variable diffusion coefficient
    const double int_a= quad_a.quad( absdet);
    if(supg_.GetSUPG())
    {   
        instat_vector_fun_ptr vel;
        if(ALE_)
        vel= Coeff::ALEVelocity;
        else
        vel= Coeff::Vel;
        QuadCL<Point3DCL> u(sit,vel,t); 
        for(int i=0; i<4; ++i)
            U_Grad[i]=dot( u, QuadCL<Point3DCL>( G[i]));
    }
    for(int i=0; i<4; ++i)
    { 
        for(int j=0; j<4; ++j)
        {
            // dot-product of the gradients

            coup[i][j]=  int_a*inner_prod( G[i], G[j])/6.0*absdet; //diffusion
            coup[i][j]+= P1DiscCL::Quad(sit, Coeff::q, i, j, 0.0)*absdet;  //reaction
            if(supg_.GetSUPG())
            {
                instat_vector_fun_ptr vel;
                if(ALE_)
                vel= Coeff::ALEVelocity;
                else
                vel= Coeff::Vel;
                QuadCL<double> res3( U_Grad[i] * U_Grad[j]);
                //SUPG stabilization
                coup[i][j]+= res3.quad(absdet)*supg_.Sta_Coeff( vel(GetBaryCenter(sit), t), Coeff::alpha(GetBaryCenter(sit), t));
            }
        }
        UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx)
                                                              : NoIdx;
    }
}

   

// mass matrix: \int_{\Omega} \alpha u  v \, dx
template<class Coeff,template <class T=double> class QuadCL>
class MassAccumulator_P1CL : public Accumulator_P1CL<Coeff,QuadCL>
{
    protected:
    typedef Accumulator_P1CL<Coeff,QuadCL> base_;
    using                           base_::MG_;
    using                           base_::BndData_; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::A_;    
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::phiQuad;
    using                           base_::coup;
    using                           base_::U_Grad;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::t;
    SUPGCL& supg_;
    bool   ALE_;
    public:
    MassAccumulator_P1CL(const MultiGridCL& MG, const BndDataCL<> * BndData, MatrixCL* Amat, VecDescCL* b, 
                   IdxDescCL& RowIdx, IdxDescCL& ColIdx, SUPGCL& supg, bool ALE, const double t)
                   :Accumulator_P1CL<Coeff,QuadCL>(MG,BndData,Amat,b,RowIdx,ColIdx,t),supg_(supg), ALE_(ALE){}
    void local_setup (const TetraCL& sit);
    void visit (const TetraCL& sit);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new MassAccumulator_P1CL ( *this); }
    
    
};


template<class Coeff,template <class T=double> class QuadCL>
void MassAccumulator_P1CL<Coeff,QuadCL>::visit (const TetraCL& sit)
{
  local_setup(sit);
  if (A_ != 0)
    base_::update_global_matrix();
  if (b_ != 0 && BndData_ != 0)      
    base_::update_coupling(sit);
}

template<class Coeff,template <class T=double> class QuadCL>
void MassAccumulator_P1CL<Coeff,QuadCL>::local_setup (const TetraCL& sit)
{

    P1DiscCL::GetGradients(G,det,sit);
    absdet= std::fabs(det);

    if(supg_.GetSUPG())
    {
        instat_vector_fun_ptr vel;
        if(ALE_)
        vel= Coeff::ALEVelocity;
        else
        vel= Coeff::Vel;
        QuadCL<Point3DCL> u(sit,vel,t);
        for(int i=0; i<4; i++)
            U_Grad[i]=dot(u, QuadCL<Point3DCL>(G[i]));
    }
    for(int i=0; i<4; ++i)
    {
      for(int j=0; j<4; ++j)
      {
        // coup[i][j]+= P1DiscCL::Quad(*sit, &Coeff::q, i, j)*absdet;
        coup[i][j]= P1DiscCL::GetMass( i, j)*absdet;
        if(supg_.GetSUPG())
        {
            instat_vector_fun_ptr vel;
            if(ALE_)
            vel= Coeff::ALEVelocity;
            else
            vel= Coeff::Vel;
            QuadCL<double> StrM(U_Grad[i]*phiQuad[j]);
            coup[i][j]+=StrM.quad(absdet)*supg_.Sta_Coeff( vel(GetBaryCenter(sit), t), Coeff::alpha(GetBaryCenter(sit), t));  //SUPG term
        }
      }
      UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx) : NoIdx;      
    }
    
}

// Convection matrix: \int_{\Omega} w \nabla u v \, dx
template<class Coeff,template <class T=double> class QuadCL>
class ConvectionAccumulator_P1CL : public Accumulator_P1CL<Coeff,QuadCL>
{
    protected:
    typedef Accumulator_P1CL<Coeff,QuadCL> base_;
    using                           base_::MG_;
    using                           base_::BndData_; 
    using                           base_::Amat_; 
    using                           base_::b_; 
    using                           base_::RowIdx_; 
    using                           base_::ColIdx_; 
    using                           base_::A_;              //Convection matrix
    using                           base_::G;
    using                           base_::det;
    using                           base_::absdet;
    using                           base_::UnknownIdx;
    using                           base_::phi;
    using                           base_::phiQuad;
    using                           base_::coup;
    using                           base_::U_Grad;
    using                           base_::lvl;
    using                           base_::idx;
    using                           base_::t;
    bool    ALE_;
    bool adjoint;
    public:
    ConvectionAccumulator_P1CL(const MultiGridCL& MG, const BndDataCL<> * BndData, MatrixCL* Amat, VecDescCL* b, 
                   IdxDescCL& RowIdx, IdxDescCL& ColIdx,  const double t_, bool ALE, bool adjoint_)
                   :Accumulator_P1CL<Coeff,QuadCL>(MG,BndData,Amat,b,RowIdx,ColIdx,t_), ALE_(ALE), adjoint(adjoint_){}
    void local_setup (const TetraCL& sit);
    void visit (const TetraCL& sit);
    virtual TetraAccumulatorCL* clone (int /*tid*/) { return new ConvectionAccumulator_P1CL ( *this); }
    
    
};


template<class Coeff,template <class T=double> class QuadCL>
void ConvectionAccumulator_P1CL<Coeff,QuadCL>::visit (const TetraCL& sit)
{
  local_setup(sit);
  if (A_ != 0)
    base_::update_global_matrix();
  if (b_ != 0 && BndData_ != 0)     
    base_::update_coupling(sit);  
}

template<class Coeff,template <class T=double> class QuadCL>
void ConvectionAccumulator_P1CL<Coeff,QuadCL>::local_setup (const TetraCL& sit)
{
    P1DiscCL::GetGradients(G,det,sit);
    absdet= std::fabs(det);

    for(int i=0; i<4; ++i)
    {
      UnknownIdx[i]= sit.GetVertex(i)->Unknowns.Exist(idx) ? sit.GetVertex(i)->Unknowns(idx) : NoIdx;
    }
    instat_vector_fun_ptr vel;
    if(ALE_)
    vel= Coeff::ALEVelocity;
    else
    vel= Coeff::Vel;
    QuadCL<Point3DCL> u(sit,vel,t);
    for(int j=0; j<4;++j)
    {
        const QuadCL<> u_Gradj( dot( u, QuadCL<Point3DCL>( G[j])));
        for(int i=0; i<4; ++i)    // assemble row i
        {
          const QuadCL<double> resq3( phiQuad[i] * u_Gradj);
          double res = resq3.quad(absdet);                  
          if (!adjoint)
            coup[i][j] = res;
          else
            coup[j][i] = res;
        }
    }
}

}//end of namespace
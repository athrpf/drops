/// \file
/// \brief Compute the interfacial tension
/// \author Hieu Nguyen, IGPM

#ifndef DROPS_SURFACETENSION_H
#define DROPS_SURFACETENSION_H

#include "num/discretize.h"
#include "num/bndData.h"

namespace DROPS
{

enum InputMethodT
/// different types of surface forces
{
    Sigma_X=0,  ///< surface tension as a function of position
    Sigma_C=1,  ///< surface tension as a function of mass concentration
    Sigma_S=2   ///< surface tension as a function of surfactant concentration
};

class SurfaceTensionCL
{
  private:
    BndDataCL<> cBnd_;
    double C_[5], cp_;                  ///< coefficients for computing surface tension with Sigma_C
    instat_scalar_fun_ptr sigma_;       ///< variable surface tension with Sigma_X
    instat_vector_fun_ptr grad_sigma_;  ///< variable surface tension gradient with Sigma_X
    InputMethodT input_;
    VecDescCL * c_;                     ///< mass concentration
    double time_;
     
    double sigma_c(double c) const;     ///< variable surface tension with Sigma_C
    double dsigma_dc(double c) const;   ///< variable surface tension gradient with Sigma_C
          
  public:
    SurfaceTensionCL ( instat_scalar_fun_ptr sigma,  instat_vector_fun_ptr grad_sigma =0)
    :  cBnd_(0),  cp_(0.), sigma_(sigma), grad_sigma_(grad_sigma), input_(Sigma_X), c_(0), time_(0.)
    { std::memset( C_, 0, 5*sizeof( double));} 
    
    SurfaceTensionCL ( double C[5], double cp=0., BndDataCL<> cBnd = BndDataCL<>( 0),  VecDescCL* c =0 , double time=0.)
    :  cBnd_(cBnd),  cp_(cp), sigma_(0), grad_sigma_(0), input_(Sigma_C), c_(c), time_(time)
    { std::memcpy( C_, C, 5*sizeof( double));}  

    void SetInputMethod(InputMethodT input) {input_=input;}
    InputMethodT GetInputMethod() {return input_;}
    void SetConcentration(VecDescCL * c) {c_=c;} 
    void SetTime(double t) {time_=t;} 
    instat_scalar_fun_ptr GetSigma() {return sigma_;}
    instat_vector_fun_ptr GetGradSigma() {return grad_sigma_;}
    void SetCoeff (double C[5], double cp) {std::memcpy( C_, C, 5*sizeof( double)); cp_= cp;}
    
    /// \brief evaluate surface tension and its gradient on triangle given by the first and second arguments 
    void ComputeSF( const TetraCL&, const BaryCoordCL * const, Quad5_2DCL<>& , Quad5_2DCL<Point3DCL>&) const; 
};

}// end of name space DROPS

#endif

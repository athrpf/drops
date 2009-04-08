//**************************************************************************
// File:    parlanczos.tpp                                                 *
// Content: lanczos algorithms                                             *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   May, 8th 2006                                                  *
//**************************************************************************
/// \author Oliver Fortmeier, SC RWTH Aachen
/// \file   parlanczos.tpp
/// \brief  Parallel lanczos algorithms for QMR itartive solver


namespace DROPS
{

/****************************************************************************
* L A N C Z O S  3  C L A S S                                               *
****************************************************************************/

template <typename Mat, typename Vec>
ParLanczos3CL<Mat,Vec>::ParLanczos3CL(const Mat& A, const Vec &v0_acc, const Vec &w0, const ExchangeCL& ex, const double BreakDownEps) :
        A_(&A), ex_(&ex),
        u_acc_(v0_acc.size()), w_acc_(w0),
        breakdownEps_(BreakDownEps),
        init_(false), breakdown_(false)
/// \param[in] A            matrix, to which Lanczos should be applied
/// \param[in] v0_acc       initial vector of basis of K(A,v0) (must be in accumulated form!)
/// \param[in] w0           initial vector of basis of K(A^T,w0) (must be in distributed form!)
/// \param[in] ex           ExchangeCL for transfer numerical data
/// \param[in] BreakDownEps if |(v,w)|<=BreakDownEps then Lanczos breaks down
{
    for (int i=0; i<2; ++i){
        v_acc_[i].resize(w0.size());
        w_[i].resize(w0.size());
    }
    w_[1]     = w0;             // remeber starting vector
    v_acc_[1] = v0_acc;
}

template <typename Mat, typename Vec>
ParLanczos3CL<Mat,Vec>::ParLanczos3CL(const double BreakDownEps) :
        breakdownEps_(BreakDownEps), init_(false), breakdown_(false)
{}

/// \brief Set nessecary information for Lanczos-Algorithm
template <typename Mat, typename Vec>
void ParLanczos3CL<Mat,Vec>::Init(const Mat& A, const Vec &v0_acc, const Vec &w0, const ExchangeCL& ex)
{
    ex_=&ex;
    A_ =&A;
    v_acc_[1]=v0_acc;
    w_[1] = w0;
    init_=false;
    breakdown_=false;
}

template <typename Mat, typename Vec>
inline void ParLanczos3CL<Mat,Vec>::SetBreakDownEps(const double eps)
{
    breakdownEps_=eps;
}

template <typename Mat, typename Vec>
inline double ParLanczos3CL<Mat,Vec>::GetBreakDownEps()
{
    return breakdownEps_;
}


/**
        One step of a modificated Lanczos-Algorithm in order to have only one sync-point.
 */
template <typename Mat, typename Vec>
bool ParLanczos3CL<Mat,Vec>::Step()
{
    if (!init_)
    {       // if not inited do a modificated step (there is no u_acc_ and no alpha_next_
        w_acc_      = w_[1];                                                    // accumulate w0
        u_acc_      = (*A_)*v_acc_[1];                                          // A*v0
        alpha_      = ex_->ParDotAcc(u_acc_,w_[1]);                             // (Av0,w0)
        ex_->Accumulate(w_acc_);                                                // send distributed entries of w0
        v_acc_.push_back( static_cast<Vec>(u_acc_-alpha_*v_acc_[1]));           // initial three term recursion for K(A,v0)
        w_.push_back( static_cast<Vec>(transp_mul((*A_),w_acc_)-alpha_*w_[1])); // initial three term recursion for K(A^T,w0)
    }
    else
    {   // here is u_acc_ the unscaled version of A*v_j ==> A*v_j = A* (\hat{v}_j/delta_j) = 1/delta_j * u_acc_
        // and alpha_next_ is ciomputed in the previous step in order to do two global reduce operations the same time
        // But w_acc_ is not accumulated yet. So do it here
        w_acc_ =w_[1];
        ex_->Accumulate(w_acc_);
        v_acc_.push_back(static_cast<Vec>((1./delta_)*u_acc_-alpha_next_*v_acc_[1]-beta_*v_acc_[0]));
        w_.push_back( static_cast<Vec>(transp_mul((*A_),w_acc_)-alpha_next_*w_[1]-delta_*w_[0]));
    }

    // calc A*v for the next step
    u_acc_ =(*A_)*v_acc_[1];
    // do two global reduce operations the same time, so collect local entries first
    val_loc_[0] = ex_->DotAcc(u_acc_, w_[1]);       // this will be alpha_next_
    val_loc_[1] = dot(v_acc_[1],w_[1]);             // = xi := (\hat{v}_{j+1},\hat{w}_{j+1})
    // now do global reduce on two doubles
    ProcCL::GlobalSum(val_loc_,val_glob_,2);

    if (!init_)         // if not inited, now the algorithm is init
        init_=true;
    else                // else store actual alpha, so the user can ask about it
        alpha_      =alpha_next_;

    delta_      = std::sqrt(std::fabs(val_glob_[1]));         // delta = std::sqrt(|xi|)
    if (delta_<=DoubleEpsC)
    {
        if (ProcCL::IamMaster())
            std::cout << "===> Breakdown of Lanczos!\n";
        breakdown_ = true;
        return false;
    }

    beta_       = val_glob_[1] / delta_;            // beta = xi/delta
    alpha_next_ = val_glob_[0] / (delta_*beta_);    // (u,\hat{w}_{j+1}) / (delta*beta)
    w_[1]      *= (1./beta_);                       // scaling of w
    v_acc_[1]  *= (1./delta_);                      // scaling of v

    return true;
}

/****************************************************************************
* L A N C Z O S  2  C L A S S                                               *
****************************************************************************/

// template <typename Mat, typename Vec>
// ParLanczos2CL<Mat,Vec>::ParLanczos2CL(const Mat& A, const Vec &v0, const Vec &w0, const ExchangeCL& ex, const double eps, const double BreakDownEps)
//  : A_(&A), ex_(&ex), v_(v0), v_acc_(v0), Av_(v0.size()), w_(w0), w_acc_(w0), p_acc_(v0.size()), Ap_(v0.size()), q_(v0.size()),
//        eps_(eps), breakdownEps_(BreakDownEps), init_(false), breakdown_(false), luckybreakdown_(false)
template <typename Mat, typename Vec, typename ExCL>
ParLanczos2CL<Mat,Vec,ExCL>::ParLanczos2CL(const Mat& A, const Vec &v0, const Vec &w0, const ExCL& ex, const double eps, const double BreakDownEps) :
    breakdownEps_(BreakDownEps)
{
    Init(A,v0,w0,ex,eps);
}

template <typename Mat, typename Vec, typename ExCL>
ParLanczos2CL<Mat,Vec,ExCL>::ParLanczos2CL(const double BreakDownEps)
    : breakdownEps_(BreakDownEps), init_(false), breakdown_(false), luckybreakdown_(false)
{}

template <typename Mat, typename Vec, typename ExCL>
void ParLanczos2CL<Mat,Vec,ExCL>::Init(const Mat& A, const Vec &v0, const Vec &w0, const ExCL& ex, const double eps)
{
    const size_t n=v0.size();
    A_=&A;
    ex_=&ex;
    v_.resize(n);     v_=v0;
    v_acc_.resize(n); v_acc_=v0;
    Av_.resize(n);    Av_=0.;
    w_.resize(n);     w_=w0;
    w_acc_.resize(n); w_acc_=w0;
    p_acc_.resize(n); p_acc_=0.;
    Ap_.resize(n);    Ap_=0.;
    q_.resize(n);     q_=0.;
    eps_=eps;
    init_=false;
    breakdown_=false;
    luckybreakdown_=false;

    val_loc_[0] = ex_->DotAcc(v_acc_,v_);
    val_loc_[1] = ex_->DotAcc(w_acc_,w_);
    Av_= (*A_)*v_acc_;
    val_loc_[2] = dot(v_,w_acc_);
    val_loc_[3] = dot(Av_,w_acc_);
    ProcCL::GlobalSum(val_loc_, val_glob_, 4);
    rho_   = val_glob_[0];
    xi_    = val_glob_[1];
    delta_ = val_glob_[2];
    sigma_ = val_glob_[3];

    if (rho_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of v!" << std::endl;
    if (xi_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of w!" << std::endl;

    rho_ = std::sqrt(std::fabs(rho_)); xi_ = std::sqrt(std::fabs(xi_));

    if (rho_<breakdownEps_)
    {breakdown_=true; luckybreakdown_=true;}
    else
    {v_ /= rho_; v_acc_ /=rho_;}

    if (xi_<breakdownEps_)
    {breakdown_=true; luckybreakdown_=true;}
    else
    {w_ /= xi_;  w_acc_ /= xi_;}
}

template <typename Mat, typename Vec, typename ExCL>
inline void ParLanczos2CL<Mat,Vec,ExCL>::SetBreakDownEps(const double eps)
{
    breakdownEps_=eps;
}

template <typename Mat, typename Vec, typename ExCL>
inline double ParLanczos2CL<Mat,Vec,ExCL>::GetBreakDownEps()
{
    return breakdownEps_;
}

template <typename Mat, typename Vec, typename ExCL>
bool ParLanczos2CL<Mat,Vec,ExCL>::Step()
{
    if (!init_){
        p_acc_ = v_acc_;
        Ap_ = Av_;
        q_ = transp_mul(*A_,w_acc_);
        q_ /=xi_;
        beta_ = sigma_/delta_;
        init_=true;
    }
    else{
        p_acc_ = v_acc_ - (delta_/(eps_*rho_)) * p_acc_;
        Ap_  = (1./rho_)*Av_ - (delta_/(eps_*rho_))*Ap_;
        q_ =  transp_mul(*A_,w_acc_) - (delta_/(eps_*xi_)) * q_;
        beta_ = sigma_/delta_ - delta_/eps_;
    }

    eps_ = (beta_*delta_)/(rho_*xi_);

    v_   = Ap_ - beta_*v_; v_acc_ = v_;
    w_   = q_ - beta_*w_;  w_acc_ = w_;

    val_loc_[0] = ex_->DotAcc(v_acc_,v_);
    Av_= (*A_)*v_acc_;
    val_loc_[1] = ex_->DotAcc(w_acc_,w_);
    val_loc_[2] = dot(v_,w_acc_);
    val_loc_[3] = dot(Av_,w_acc_);
    ProcCL::GlobalSum(val_loc_, val_glob_, 4);

    rho_   = val_glob_[0];
    xi_    = val_glob_[1];
    delta_ = val_glob_[2];
    sigma_ = val_glob_[3];

    if (rho_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of v!" << std::endl;
    if (xi_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of w!" << std::endl;

    rho_ = std::sqrt(std::fabs(rho_));
    xi_ = std::sqrt(std::fabs(xi_));

    bool ret = true;

    bool check_delta = std::fabs(delta_) < breakdownEps_,
         check_rho   =rho_<breakdownEps_,
         check_xi    =xi_<breakdownEps_;

    breakdown_ = check_delta || check_rho || check_xi;

    if (breakdown_){
        luckybreakdown_ = check_rho || check_xi;
        if (ProcCL::IamMaster())
            std::cout << "===> "
                    << (luckybreakdown_? "lucky" : "serious")
                    << " Breakdown of coupled Lanczos, because "
                    << (check_rho ? "v is nearly zero" : "")
                    << (check_xi  ? "w is nearly zero" : "")
                    << ", delta " <<delta_<< ", rho " <<rho_<<", xi "<<xi_
                    << std::endl;

        ret = false;
    }

    if (!check_rho){
        v_ /= rho_; v_acc_ /=rho_;
    }

    if (!check_xi){
        w_ /= xi_;  w_acc_ /= xi_;
    }

    return ret;
}


template<typename Mat, typename Vec, typename PC, typename ExCL>
ParPreLanczos2CL<Mat,Vec,PC,ExCL>::ParPreLanczos2CL(const Mat &A, const Vec &v0, const Vec &w0,
                                                    const PC &pc, const ExCL &ex,
                                                    const double eps, const double BreakDownEps) :
        pc_(pc), breakdownEps_(BreakDownEps)
{
    Init(A,v0,w0,ex,eps);
}

template <typename Mat, typename Vec, typename PC, typename ExCL>
ParPreLanczos2CL<Mat,Vec,PC,ExCL>::ParPreLanczos2CL(const PC &pc, const double BreakDownEps)
    : pc_(pc), breakdownEps_(BreakDownEps), init_(false), breakdown_(false), luckybreakdown_(false)
{}

template <typename Mat, typename Vec, typename PC, typename ExCL>
void ParPreLanczos2CL<Mat,Vec,PC,ExCL>::Init(const Mat& A, const Vec &v0, const Vec &w0, const ExCL& ex, const double eps)
{
    const size_t n=v0.size();
    A_=&A;
    ex_=&ex;
    v_.resize(n);     v_=v0;
    v_acc_.resize(n); v_acc_=v0;
    Av_.resize(n);    Av_=0.;
    w_.resize(n);     w_=w0;
    w_acc_.resize(n); w_acc_=w0;
    p_acc_.resize(n); p_acc_=0.;
    Ap_.resize(n);    Ap_=0.;
    q_.resize(n);     q_=0.;
    eps_=eps;
    init_=false;
    breakdown_=false;
    luckybreakdown_=false;

    val_loc_[0] = ex_->DotAcc(v_acc_,v_);
    val_loc_[1] = ex_->DotAcc(w_acc_,w_);
    pc_.Apply((*A_), Av_ ,(*A_)*v_acc_);
    val_loc_[2] = dot(v_,w_acc_);
    val_loc_[3] = dot(Av_,w_acc_);
    ProcCL::GlobalSum(val_loc_, val_glob_, 4);
    rho_   = val_glob_[0];
    xi_    = val_glob_[1];
    delta_ = val_glob_[2];
    sigma_ = val_glob_[3];

    if (rho_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of v!" << std::endl;
    if (xi_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of w!" << std::endl;

    rho_ = std::sqrt(std::fabs(rho_)); xi_ = std::sqrt(std::fabs(xi_));

    if (rho_<breakdownEps_)
    {breakdown_=true; luckybreakdown_=true;}
    else
    {v_ /= rho_; v_acc_ /=rho_;}

    if (xi_<breakdownEps_)
    {breakdown_=true; luckybreakdown_=true;}
    else
    {w_ /= xi_;  w_acc_ /= xi_;}
}

template <typename Mat, typename Vec, typename PC, typename ExCL>
inline void ParPreLanczos2CL<Mat,Vec,PC,ExCL>::SetBreakDownEps(const double eps)
{
    breakdownEps_=eps;
}

template <typename Mat, typename Vec, typename PC, typename ExCL>
inline double ParPreLanczos2CL<Mat,Vec,PC,ExCL>::GetBreakDownEps()
{
    return breakdownEps_;
}

template <typename Mat, typename Vec, typename PC, typename ExCL>
bool ParPreLanczos2CL<Mat,Vec,PC,ExCL>::Step()
{
    if (!init_){
        p_acc_ = v_acc_;
        Ap_ = Av_;
        VectorCL tmp(w_acc_.size());
        pc_.transp_Apply(*A_,tmp, w_acc_);
        q_ = transp_mul(*A_,tmp);
        q_ /=xi_;
        beta_ = sigma_/delta_;
        init_=true;
    }
    else{
        p_acc_ = v_acc_ - (delta_/(eps_*rho_)) * p_acc_;
        Ap_  = (1./rho_)*Av_ - (delta_/(eps_*rho_))*Ap_;

        VectorCL tmp(w_acc_.size());
        pc_.transp_Apply(*A_,tmp, w_acc_);
        tmp = transp_mul(*A_,tmp);

        q_ =  tmp - (delta_/(eps_*xi_)) * q_;
        beta_ = sigma_/delta_ - delta_/eps_;
    }

    eps_ = (beta_*delta_)/(rho_*xi_);

    v_   = Ap_ - beta_*v_; v_acc_ = v_;
    w_   = q_ - beta_*w_;  w_acc_ = w_;

    val_loc_[0] = ex_->DotAcc(v_acc_,v_);
    pc_.Apply(*A_,Av_,(*A_)*v_acc_);
    val_loc_[1] = ex_->DotAcc(w_acc_,w_);
    val_loc_[2] = dot(v_,w_acc_);
    val_loc_[3] = dot(Av_,w_acc_);
    ProcCL::GlobalSum(val_loc_, val_glob_, 4);

    rho_   = val_glob_[0];
    xi_    = val_glob_[1];
    delta_ = val_glob_[2];
    sigma_ = val_glob_[3];

    if (rho_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of v!" << std::endl;
    if (xi_<0)
        std::cout << "ParLanczos2CL::Step(): negative norm of w!" << std::endl;

    rho_ = std::sqrt(std::fabs(rho_));
    xi_ = std::sqrt(std::fabs(xi_));

    bool ret = true;

    bool check_delta = std::fabs(delta_) < breakdownEps_,
         check_rho   = rho_<breakdownEps_,
         check_xi    = xi_<breakdownEps_;

    breakdown_ = check_delta || check_rho || check_xi;

    if (breakdown_){
        luckybreakdown_ = check_rho || check_xi;
        if (ProcCL::IamMaster())
            std::cout << "===> "
                    << (luckybreakdown_? "lucky" : "serious")
                    << " Breakdown of coupled Lanczos, because "
                    << (check_delta ? "(v,w) is nearly zero" : "")
                    << (check_rho ? "v is nearly zero" : "")
                    << (check_xi  ? "w is nearly zero" : "")
                    << ", delta " <<delta_<< ", rho " <<rho_<<", xi "<<xi_
                    << std::endl;

        ret = false;
    }

    if (!check_rho){
        v_ /= rho_; v_acc_ /=rho_;
    }

    if (!check_xi){
        w_ /= xi_;  w_acc_ /= xi_;
    }

    return ret;
}

} // end of namespace DROPS

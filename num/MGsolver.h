/// \file
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Maxim Larin, Patrick Esser, IGPM

#ifndef DROPS_MGSOLVER_H
#define DROPS_MGSOLVER_H

#include "misc/problem.h"
#include "num/solver.h"

#include <list>


namespace DROPS
{

struct MGLevelDataCL  // data for one triang level
{
    IdxDescCL Idx, IdxPr;          // index description for velocity and pressure
    MatDescCL   A,     B,     BT;  // saddle-point matrix
    MatDescCL   Mpr, Mvel;         // mass matrices
    MatDescCL   P,   PPr;          // prolongation matrices for velocity and pressure
    MatDescCL   AN;                // A + alpha*N
    MatrixCL*   ABlock;
};

class MGDataCL : public std::list<MGLevelDataCL>
{
  private:
    bool StokesMG_;

  public:
    MGDataCL(int n=0) : std::list<MGLevelDataCL>(n), StokesMG_(false) {}
    void RemoveCoarseResetFinest() {
        if (this->empty()) return;
        //RemoveCoarse
        MGDataCL::iterator it=this->end();
        --it;
        this->erase(this->begin(), it);
        //ResetFinest
        this->begin()->A.Data.clear();
        this->begin()->P.Data.clear();
        this->begin()->B.Data.clear();
        this->begin()->BT.Data.clear();
        this->begin()->Mpr.Data.clear();
        this->begin()->Mvel.Data.clear();
        this->begin()->PPr.Data.clear();
        this->begin()->AN.Data.clear();
        this->begin()->ABlock = &this->begin()->A.Data;
    }
    bool StokesMG() {return StokesMG_;}
    void SetStokesMG(bool full) {StokesMG_=full;}
};

typedef MGDataCL::iterator       MGDataIterCL;
typedef MGDataCL::const_iterator const_MGDataIterCL;

// Multigrid method, V-cycle, beginning from level 'fine'
// numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
//      after 'numLevel' visited levels or if #Unknowns <= 'numUnknDirect'
// If one of the parameters is -1, it will be neglected.
// If the coarsest level 'begin' has been reached, the direct solver is used too.
// NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)
template<class SmootherCL, class DirectSolverCL>
void MGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b,
          const SmootherCL& Smoother, Uint smoothSteps,
          DirectSolverCL& Solver, int numLevel, int numUnknDirect);


// checks  A_coarse= PT * A_fine * P on each level
void CheckMGData( const_MGDataIterCL begin, const_MGDataIterCL end);


// Uses MGM for solving to tolerance tol or until maxiter iterations are reached.
// The error is measured as two-norm of dx for residerr=false, of Ax-b for residerr=true.
// sm controls the number of smoothing steps, lvl the number of used levels
template<class SmootherCL, class DirectSolverCL>
void MG(const MGDataCL& MGData, const SmootherCL&, DirectSolverCL&, VectorCL& x, const VectorCL& b,
        int& maxiter, double& tol, const bool residerr= true, Uint sm=1, int lvl=-1);


template<class SmootherT, class DirectSolverT>
class MGSolverCL : public SolverBaseCL
{
  private:
    const MGDataCL&      mgdata_;
    const SmootherT&     smoother_;
    DirectSolverT&       directSolver_;
    const bool residerr_;
    Uint smoothSteps_;
    int usedLevels_;

  public:
    MGSolverCL( const MGDataCL& mgdata, const SmootherT& sm, DirectSolverT& ds, int maxiter, double tol,
               const bool residerr= true, Uint smsteps= 1, int lvl= -1 )
        : SolverBaseCL(maxiter,tol), mgdata_(mgdata), smoother_(sm), directSolver_(ds),
          residerr_(residerr), smoothSteps_(smsteps), usedLevels_(lvl) {}

    void Solve(const MatrixCL& /*A*/, VectorCL& x, const VectorCL& b)
    {
        _res=  _tol;
        _iter= _maxiter;
        MG( mgdata_, smoother_, directSolver_, x, b, _iter, _res, residerr_, smoothSteps_, usedLevels_);
    }
};

// Multigrid method, V-cycle, beginning from level 'fine' 
// numLevel and numUnknDirect specify, when the direct solver 'Solver' is used:
//      after 'numLevel' visited levels or if #Unknowns <= 'numUnknDirect'
// If one of the parameters is -1, it will be neglected. 
// If the coarsest level 'begin' has been reached, the direct solver is used too.
// NOTE: Assumes, that the levels are stored in an ascending order (first=coarsest, last=finest)
template<class StokesSmootherCL, class StokesDirectSolverCL>
void StokesMGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& u, VectorCL& p,
                const VectorCL& b, const VectorCL& c,  const StokesSmootherCL& Smoother, Uint smoothSteps,
                Uint cycleSteps, StokesDirectSolverCL& Solver, int numLevel, int numUnknDirect);

// checks Stokes-MG-Data
void CheckStokesMGData( const_MGDataIterCL begin, const_MGDataIterCL end);

//===================================
// definition of template functions
//===================================

template<class SmootherCL, class DirectSolverCL>
void
MGM(const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& x, const VectorCL& b,
    const SmootherCL& Smoother, Uint smoothSteps,
    DirectSolverCL& Solver, int numLevel, int numUnknDirect)
// Multigrid method, V-cycle. If numLevel==0 or #Unknowns <= numUnknDirect,
// the direct solver Solver is used.
// If one of the parameters is -1, it will be neglected.
// If MGData.begin() has been reached, the direct solver is used too.
{
    const_MGDataIterCL coarse= fine;
    --coarse;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) )
       || fine==begin)
    { // use direct solver
        // NOTE: Use relative residual-measurement, as otherwise the accuracy and
        // correctness depend on the scaling of the matrix and geometry.
        // This has bitten us in e.g. levelset/mzelle_instat.cpp.
        Solver.Solve( *fine->ABlock, x, b);
        std::cerr << "MGM: direct solver: iterations: " << Solver.GetIter()
                  << "\tresiduum: " << Solver.GetResid() << '\n';
        return;
    }
    VectorCL d(coarse->Idx.NumUnknowns), e(coarse->Idx.NumUnknowns);
    // presmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine->ABlock, x, b);
    // restriction of defect
    d= transp_mul( fine->P.Data, VectorCL( b - *fine->ABlock*x));
    // calculate coarse grid correction
    MGM( begin, coarse, e, d, Smoother, smoothSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    // add coarse grid correction
    x+= fine->P.Data * e;
    // postsmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine->ABlock, x, b);
}


template<class SmootherCL, class DirectSolverCL>
void MG(const MGDataCL& MGData, const SmootherCL& smoother, DirectSolverCL& solver, 
        VectorCL& x, const VectorCL& b,
        int& maxiter, double& tol, const bool residerr, Uint sm, int lvl)
{
    const_MGDataIterCL finest= --MGData.end();
    double resid= -1;
    double old_resid;
    VectorCL tmp;
    if (residerr == true) {
        resid= norm( b - *finest->ABlock * x);
        //std::cerr << "initial residual: " << resid << '\n';
    }
    else
        tmp.resize( x.size());

    int it;
    for (it= 0; it<maxiter; ++it) {
        if (residerr == true) {
            if (resid <= tol) break;
        }
        else tmp= x;
        MGM( MGData.begin(), finest, x, b, smoother, sm, solver, lvl, -1);
        if (residerr == true) {
            old_resid= resid;
            resid= norm( b - *finest->ABlock * x);
//            std::cerr << "iteration: " << it  << "\tresidual: " << resid;
//            std::cerr << "\treduction: " << resid/old_resid;
//            std::cerr << '\n';
        }
        else if ((resid= norm( tmp - x)) <= tol) break;
    }
    maxiter= it;
    tol= resid;
}

template<class StokesSmootherCL, class StokesDirectSolverCL>
void StokesMGM( const const_MGDataIterCL& begin, const const_MGDataIterCL& fine, VectorCL& u, VectorCL& p,
                const VectorCL& b, const VectorCL& c,  const StokesSmootherCL& Smoother, Uint smoothSteps,
                Uint cycleSteps, StokesDirectSolverCL& Solver, int numLevel, int numUnknDirect)
//
// Multigrid method for Stokes problem: V-cycle. 
// If numLevel==0 or #Unknowns <= numUnknDirect, the direct solver is used.
// If one of the parameters is -1, it will be neglected.
// If MGData.begin() has been reached, the direct solver is used too.
//
{
    const_MGDataIterCL coarse= fine;
    --coarse;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : u.size() <= static_cast<Uint>(numUnknDirect) ) 
       || fine==begin)
    {
        // use direct solver
        std::cout << "P2P1:StokesMGM: use direct solver " << std::endl;
        Solver.Solve( *fine->ABlock, fine->B.Data, u, p, b, c);
        std::cerr << "P2P1:StokesMGM: direct solver: iter: " << Solver.GetIter() << "\tresid: " << Solver.GetResid() << std::endl;
        return;
    }

    VectorCL du(coarse->Idx.NumUnknowns), dp(coarse->IdxPr.NumUnknowns);
    VectorCL eu(coarse->Idx.NumUnknowns), ep(coarse->IdxPr.NumUnknowns);
    // presmoothing
//    std::cout << "P2P1:StokesMGM: presmoothing " << smoothSteps << " steps " << std::endl;
    for (Uint i=0; i<smoothSteps; ++i) {
        Smoother.Apply( *fine->ABlock, fine->B.Data, fine->BT.Data, fine->Mpr.Data, u, p, b, c );
    }
    // restriction of defect
//    std::cout << "P2P1:StokesMGM: restriction of defect " << std::endl;
    du = transp_mul( fine->P.Data, VectorCL(b - (*fine->ABlock*u + transp_mul( fine->B.Data, p ) )) );
    dp = transp_mul( fine->PPr.Data,  VectorCL(c - fine->B.Data*u ));

    // calculate coarse grid correction
//    std::cout << "P2P1:StokesMGM: calculate coarse grid correction    " << cycleSteps << " times " << std::endl;
    eu=0;
    ep=0;
    for (Uint i=0; i<cycleSteps; ++i) {
      StokesMGM( begin, coarse, eu, ep, du, dp, Smoother, smoothSteps, cycleSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    }
    // add coarse grid correction
//    std::cout << "P2P1:StokesMGM: add coarse grid correction " << std::endl;
    u+= fine->P.Data * eu;
    p+= fine->PPr.Data  * ep;
    // postsmoothing
//    std::cout << "P2P1:StokesMGM: postsmoothing " << smoothSteps << " steps " << std::endl;
    for (Uint i=0; i<smoothSteps; ++i) {
        Smoother.Apply( *fine->ABlock, fine->B.Data, fine->BT.Data, fine->Mpr.Data, u, p, b, c );
    }
}

class PVankaSmootherCL
{
  private:
    int vanka_method_;
    double tau_;
// choose the Vanka smoother:
//                            0 - diagonal
//                            2 - full (LR)
//                            3 - Gauss-Seidel
//                            4 - SSOR (symmetrical Gauss-Seidel)

    template <typename Mat>
    DMatrixCL<double>  SetupLocalProblem (Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B) const;
    template <typename Mat, typename Vec>
    void DiagSmoother(Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;
    template <typename Mat, typename Vec>
    void GSSmoother(Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;
    template <typename Mat, typename Vec>
    void LRSmoother(Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B, Vec& v) const;

  public:
    PVankaSmootherCL(int vanka_method=0, double tau=1.0) : vanka_method_(vanka_method), tau_(tau) {}
    template <typename Mat, typename Vec>
    void Apply( const Mat& A, const Mat& B, const Mat& BT, const Mat&,
                Vec& u, Vec& p, const Vec& f, const Vec& g ) const;
    void SetVankaMethod( int method) {vanka_method_ = method;}
    int  GetVankaMethod()            {return vanka_method_;}
};


// constructs local saddle point matrix
// TODO effizenter implementieren (CRS-Struktur ausnutzen)
template <typename Mat>
DMatrixCL<double> PVankaSmootherCL::SetupLocalProblem (Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B) const
{
    const size_t dim = NodeListVel.size();
    DMatrixCL<double> M(dim+1, dim+1);
    for ( size_t i=0; i<dim; ++i )
    {
        const int irow = NodeListVel[i];
        for ( size_t j=0; j<dim; ++j )
        {
            const int jcol = NodeListVel[j];
            M(i,j) = A(irow,jcol);
        }
        M(i,dim) = M(dim,i) = B(id,irow);
    }
    M(dim,dim) = 0.0;
    return M;
}

template <typename Mat, typename Vec>
void PVankaSmootherCL::DiagSmoother(Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B, Vec& v) const
{
    const size_t dim = v.size()-1;
    Vec DiagA(dim);
    Vec Block(dim);
    // Direct solver: constructing the local problem Mx=v
    for ( size_t i=0; i<dim; ++i )
    {
        const int irow = NodeListVel[i];
        DiagA[i] = A(irow,irow);
    }

    const size_t n= B.row_beg( id + 1) - B.row_beg( id);
    std::memcpy( Addr(Block), B.raw_val() + B.row_beg( id), n*sizeof(typename Vec::value_type));
    // Solve Mx=v

    double S=0.0;
    for ( size_t i=0; i<dim; ++i ){
        v[i] /= DiagA[i];
        v[dim] -= Block[i]*v[i];
        S += Block[i]*Block[i]/DiagA[i];
    }
    v[dim] /= -S;
    for ( size_t i=0; i<dim; ++i )
        v[i]-=v[dim]*Block[i]/DiagA[i];
}

template <typename Mat, typename Vec>
void GaussSeidel(const Mat& M, Vec& x, const Vec& rhs)
{
    const size_t dim = rhs.size();
    x = rhs;
    for ( size_t i=0; i<dim; ++i )
    {
        for ( size_t j=0; j<i; ++j )
            x[i]-=M(i,j)*x[j];
        x[i]/=M(i,i);
    }
}

template <typename Mat, typename Vec>
void SymmetricGaussSeidel(const Mat& M, Vec& x, const Vec& rhs)
{
    const size_t dim = rhs.size();
    x = rhs;
    for ( size_t i=0; i<dim; ++i )
    {
        for ( size_t j=0; j<i; ++j )
            x[i]-=M(i,j)*x[j];
        x[i]/=M(i,i);
    }

    double tmp=0.;
    for ( size_t i=dim-1; i<dim; --i)
    {
        tmp=rhs[i];
        for (size_t j=0; j<dim; ++j)
            tmp-= M(i,j)*x[j];
        x[i]+=tmp/M(i,i);
    }
}

template <typename Mat, typename Vec>
void Jacobi(const Mat& M, Vec& x, const Vec& rhs)
{
    const size_t dim = rhs.size();
    x = rhs;
    for ( size_t i=0; i<dim; ++i )
        x[i]/=M(i,i);
}


// exact schur method with A^{-1} = one step (symmetric) Gauss Seidel
template <typename Mat, typename Vec>
void PVankaSmootherCL::GSSmoother(Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B, Vec& v) const
{
    const int dim = v.size() - 1;
    DMatrixCL<double> M(SetupLocalProblem(id, NodeListVel, A, B));
    Vec rhs (dim), x(dim), y(dim);

    for ( int i=0; i<dim; ++i )
        rhs[i]=M(i,dim);

    //A*y=B^T
    //Jacobi(M, y, rhs);
    switch (vanka_method_)
    {
        case 3 : GaussSeidel(M, y, rhs);          break;
        case 4 : SymmetricGaussSeidel(M, y, rhs); break;
    }

    double schur=0.0;
    for ( int i=0; i<dim; ++i )
        schur += M(dim,i)*y[i];

    for ( int i=0; i<dim; ++i )
        rhs[i]=v[i];

    //Jacobi(M, x, rhs);
    switch (vanka_method_)
    {
        case 3 : GaussSeidel(M, x, rhs);          break;
        case 4 : SymmetricGaussSeidel(M, x, rhs); break;
    }

    double rhsd=0.0;
    for ( int i=0; i<dim; ++i )
        rhsd += M(dim,i)*x[i];
    rhsd -= v[dim];
    v[dim] = rhsd/schur;

    for ( int i=0; i<dim; ++i )
        v[i]=x[i]-v[dim]*y[i];//M(i,dim);

}

template <typename Mat, typename Vec>
void PVankaSmootherCL::LRSmoother(Uint id, std::vector<size_t>& NodeListVel, const Mat& A, const Mat& B, Vec& v) const
{
    DMatrixCL<double> M(SetupLocalProblem(id, NodeListVel, A, B));
    gauss_pivot (M, v);
}

template <typename Mat, typename Vec>
void
PVankaSmootherCL::Apply(const Mat& A, const Mat& B, const Mat& BT, const Mat&, Vec& x, Vec& y, const Vec& f, const Vec& g) const
{
    std::vector<size_t> NodeListVel;
    NodeListVel.reserve(200);
    size_t dim;
    // Cycle by pressure variables
    for ( size_t id=0; id<g.size(); ++id )
    {
        dim = B.row_beg(id+1) - B.row_beg(id);
        Vec v(dim+1);

        NodeListVel.resize(dim);

        // copy column indices of nonzeros from row id of B to NodeListVel
        std::memcpy( &NodeListVel[0], B.raw_col() + B.row_beg( id), dim*sizeof(size_t));

        // v = resid of (A & B^T \\ B & 0) * (x \\ y) - (f \\ g) for the rows NodeListVel, id
        for ( size_t i=0; i<dim; ++i )
        {
            const int irow = NodeListVel[i];
            double sum=0;
            for ( Uint k=A.row_beg( irow ); k<A.row_beg( irow+1 ); ++k )
                sum+= A.val(k)*x[A.col_ind(k)];

            for ( Uint k=BT.row_beg( irow ); k<BT.row_beg( irow+1 ); ++k )
                sum+= BT.val(k)*y[BT.col_ind(k)];
            sum-=f[irow];
            v[i]=sum;
        }

        double sump=0;
        for ( size_t k=B.row_beg( id ); k<B.row_beg( id+1 ); ++k )
            sump+= B.val(k)*x[B.col_ind(k)];
        sump-=g[id];
        v[dim]=sump;

        // smoothing
        switch (vanka_method_) {
            case 0 : {
                DiagSmoother(id, NodeListVel, A, B, v);
            } break;
            case 2: {
                LRSmoother  (id, NodeListVel, A, B, v);
            } break;
            case 3: case 4 : {
                GSSmoother  (id, NodeListVel, A, B, v);
            } break;
        }

        // Cycle by local unknowns: correction of the approximation x
        for ( size_t i=0; i<dim; ++i )
        {
            const size_t indvel = NodeListVel[i];
            x[indvel] -= tau_ * v[i];
        }
        y[id] -= tau_ * v[dim];
    }
}

class BSSmootherCL
{
  private:
    int    maxit_;
    double red_, omega_;

  public:
    BSSmootherCL( int maxit = 20, double red = 2e-1, double omega = 2.0) : maxit_(maxit), red_(red), omega_(omega) {};
    template <typename Mat, typename Vec>
    void Apply( const Mat& A, const Mat& B, const Mat&, const Mat& M,
                Vec& u, Vec& p, const Vec& f, const Vec& g ) const;
};

template <typename Mat, typename Vec>
void BSSmootherCL::Apply( const Mat& A, const Mat& B, const Mat&, const Mat& M, Vec& u, Vec& p, const Vec& f, const Vec& g ) const
{
     const size_t n= f.size();
     const size_t m= g.size();

     double actualtol = 1;

// compute the residual of the Stokes problem

     Vec ru(A * u + transp_mul(B, p) - f);
     Vec rp(B * u                    - g);

     Vec x(n), y(m), z(n);
     Vec DiagA(A.GetDiag()*omega_);

// forward substitution

     x = ru/DiagA;

     y = rp - B * x;

     Vec xk(m);
     Vec pk(m), zk(m), qk(m);

     Vec r (-y);
     double rho, rho_1= 0.0, resid= norm_sq(r);
     double rhsnorm = resid;

    for (int i=1; i<=maxit_; ++i) // CG for schur complement
    {
        zk=r;
        rho= dot(r,zk);
        if (i == 1)
            pk= zk;
        else
            z_xpay(pk, zk, (rho/rho_1), pk); // p= z + (rho/rho_1)*p;
        z=transp_mul(B, pk);
        z/=DiagA;
        qk = B * z;
        const double alpha= rho/dot(pk,qk);
        axpy(alpha, pk, xk);                // x+= alpha*p;
        axpy(-alpha, qk, r);                // r-= alpha*q;

        resid= norm_sq(r);
        actualtol= resid/rhsnorm;
        if (actualtol<=red_)
            break;
        rho_1= rho;
    }
// backward substitution

    z = transp_mul(B, xk);
    z /= DiagA;

    u -= x-z;
    p -= xk;
// p must be in L_2^0(\Omega) TODO: other boundary conditions
    Vec ones( 1.0, p.size() );
    Vec oneM= M * ones;
    double coef= dot( oneM, p ) / dot( oneM, ones );
    p -= coef*ones;
}


} // end of namespace DROPS

#endif

/// \file MGsolver.tpp
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author LNM RWTH Aachen: Sven Gross,  Patrick Esser, Joerg Grande, Maxim Larin, Volker Reichelt; SC RWTH Aachen:

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

namespace DROPS {

template <class SmootherCL, class DirectSolverCL, class ProlongationIteratorT>
void
MGM(const MLMatrixCL::const_iterator& begin, const MLMatrixCL::const_iterator& fine,
     const ProlongationIteratorT& P, VectorCL& x, const VectorCL& b,
     const SmootherCL& Smoother, Uint smoothSteps,
     DirectSolverCL& Solver, int numLevel, int numUnknDirect)
{
    MLMatrixCL::const_iterator coarse= fine;
    ProlongationIteratorT      coarseP= P;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) )
       || fine==begin)
    { // use direct solver
        Solver.Solve( *fine, x, b);
/*        std::cout << "MGM: direct solver: iterations: " << Solver.GetIter()
                  << "\tresiduum: " << Solver.GetResid() << '\n';*/
        return;
    }
    --coarse;
    --coarseP;
    VectorCL d( coarse->num_cols()), e( coarse->num_cols());
    // presmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine, x, b);
    // restriction of defect
    d= transp_mul( *P, VectorCL( b - *fine*x));
    // calculate coarse grid correction
    MGM( begin, coarse, coarseP, e, d, Smoother, smoothSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    // add coarse grid correction
    x+= (*P) * e;
    // postsmoothing
    for (Uint i=0; i<smoothSteps; ++i) Smoother.Apply( *fine, x, b);
}

template<class SmootherCL, class DirectSolverCL, class ProlongationT>
void MG(const MLMatrixCL& MGData, const ProlongationT& Prolong, const SmootherCL& smoother,
        DirectSolverCL& solver, VectorCL& x, const VectorCL& b, int& maxiter, double& tol,
        const bool residerr, Uint sm, int lvl)
{
    MLMatrixCL::const_iterator finest= MGData.GetFinestIter();
    typename ProlongationT::const_iterator finestProlong= Prolong.GetFinestIter();
    double resid= -1;
//    double old_resid;
    VectorCL tmp;
    if (residerr == true) {
        resid= norm( b - *finest * x);
        //std::cout << "initial residual: " << resid << '\n';
    }
    else
        tmp.resize( x.size());

    int it;
    for (it= 0; it<maxiter; ++it) {
        if (residerr == true) {
            if (resid <= tol) break;
        }
        else tmp= x;
        MGM( MGData.begin(), finest, finestProlong, x, b, smoother, sm, solver, lvl, -1);
        if (residerr == true) {
//            old_resid= resid;
            resid= norm( b - *finest * x);
//            std::cout << "iteration: " << it  << "\tresidual: " << resid;
//            std::cout << "\treduction: " << resid/old_resid;
//            std::cout << '\n';
        }
        else if ((resid= norm( tmp - x)) <= tol) break;
    }
    maxiter= it;
    tol= resid;
}

template<class StokesSmootherCL, class StokesDirectSolverCL, class ProlongItT1, class ProlongItT2>
void StokesMGM( const MLMatrixCL::const_iterator& beginA,  const MLMatrixCL::const_iterator& fineA,
                const MLMatrixCL::const_iterator& fineB,   const MLMatrixCL::const_iterator& fineBT, 
                const MLMatrixCL::const_iterator& fineprM, const ProlongItT1& PVel,
                const ProlongItT2& PPr, VectorCL& u, VectorCL& p, const VectorCL& b,
                const VectorCL& c, const StokesSmootherCL& Smoother, Uint smoothSteps, Uint cycleSteps,
                StokesDirectSolverCL& Solver, int numLevel, int numUnknDirect)
{
    MLMatrixCL::const_iterator coarseA   = fineA;
    MLMatrixCL::const_iterator coarseB   = fineB;
    MLMatrixCL::const_iterator coarseBT  = fineBT;
    MLMatrixCL::const_iterator coarseprM = fineprM;
    ProlongItT1 coarsePVel= PVel;
    ProlongItT2 coarsePPr = PPr;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : u.size() <= static_cast<Uint>(numUnknDirect) )
       || fineA==beginA)
    {
        // use direct solver
        std::cout << "P2P1:StokesMGM: use direct solver " << std::endl;
        Solver.Solve( *fineA, *fineB, u, p, b, c);
        //std::cout << "P2P1:StokesMGM: direct solver: iter: " << Solver.GetIter() << "\tresid: " << Solver.GetResid() << std::endl;
        return;
    }
    --coarseA;
    --coarseB;
    --coarseBT;
    --coarseprM;
    --coarsePVel;
    --coarsePPr;

    // presmoothing
//    std::cout << "P2P1:StokesMGM: presmoothing " << smoothSteps << " steps " << std::endl;
    for (Uint i=0; i<smoothSteps; ++i) {
        Smoother.Apply( *fineA, *fineB, *fineBT, *fineprM, u, p, b, c );
    }
    // restriction of defect
    //std::cout << "P2P1:StokesMGM: restriction of defect " << std::endl;

    VectorCL du (transp_mul( *PVel, VectorCL(b - (*fineA * u + transp_mul( *fineB, p ) )) ));
    VectorCL dp (transp_mul( *PPr,  VectorCL(c - *fineB * u )));

    // calculate coarse grid correction
//    std::cout << "P2P1:StokesMGM: calculate coarse grid correction    " << cycleSteps << " times " << std::endl;
    VectorCL eu ( du.size());
    VectorCL ep ( dp.size());
    for (Uint i=0; i<cycleSteps; ++i) {
      StokesMGM( beginA, coarseA, coarseB, coarseBT, coarseprM, coarsePVel, coarsePPr, eu, ep, du, dp, Smoother,
                 smoothSteps, cycleSteps, Solver, (numLevel==-1 ? -1 : numLevel-1), numUnknDirect);
    }
    // add coarse grid correction
//    std::cout << "P2P1:StokesMGM: add coarse grid correction " << std::endl;
    u+= *PVel * eu;
    p+= *PPr * ep;
    // postsmoothing
//    std::cout << "P2P1:StokesMGM: postsmoothing " << smoothSteps << " steps " << std::endl;
    for (Uint i=0; i<smoothSteps; ++i) {
        Smoother.Apply( *fineA, *fineB, *fineBT, *fineprM, u, p, b, c );
    }
}

// TODO effizenter implementieren (CRS-Struktur ausnutzen)
template <typename Mat>
DMatrixCL<double> PVankaSmootherCL::SetupLocalProblem (Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B) const
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

inline void
Solve2x2 (const SMatrixCL<2,2>& A, double& x0, double& x1)
{
    const double det( A(0,0)*A(1,1) - A(1,0)*A(0,1));
    const double b0( x0);
    x0= (A(1,1)*x0-A(0,1)*x1)/det;
    x1= (A(0,0)*x1-A(1,0)*b0)/det;
}

template <typename Mat, typename Vec>
void PVankaSmootherCL::DiagSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v, size_t id2) const
{
    const size_t dim = v.size() - (id2 == NoIdx ? 1 : 2);
    Vec DiagA( dim);
    Vec Block( dim);
    Vec Block2( dim);

   // Direct solver: constructing the local problem Mx=v
    for (size_t i= 0; i < dim; ++i) {
        const int irow = NodeListVel[i];
        DiagA[i] = A( irow, irow);
    }

    const size_t n= B.row_beg( id + 1) - B.row_beg( id);
    std::memcpy( Addr(Block), B.raw_val() + B.row_beg( id), n*sizeof(typename Vec::value_type));
    if (id2 != NoIdx) { // copy the row for id2
        const size_t* idbeg=  B.GetFirstCol( id);
        const size_t* id2beg= B.GetFirstCol( id2);
        const double* idvalbeg=  B.GetFirstVal( id);
        const double* id2valbeg= B.GetFirstVal( id2);
        for (size_t pos= 0; idbeg < B.GetFirstCol( id + 1); ++pos, ++idbeg, ++idvalbeg) {
            if (*idbeg < *id2beg)
                Block2[pos]= 0.;
            else {
                Block2[pos]= *id2valbeg++;
                ++id2beg;
            }
        }
    }
    
    // Solve Mx=v
    if (id2 == NoIdx) {
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
    else {
        SMatrixCL<2,2> S;
        for ( size_t i=0; i<dim; ++i ){
            S( 0, 0) += Block[i]*Block[i]/DiagA[i];
            S( 1, 0) += Block[i]*Block2[i]/DiagA[i];
            S( 1, 1) += Block2[i]*Block2[i]/DiagA[i];
        }
        S( 0, 1)= S( 1, 0);
        S*= -1.;
        v[std::slice( 0, dim, 1)]/= DiagA;
        v[dim]    -= dot( Block,  Vec( v[std::slice( 0, dim, 1)]));
        v[dim + 1]-= dot( Block2, Vec( v[std::slice( 0, dim, 1)]));
        Solve2x2( S, v[dim], v[dim + 1]);
        for (size_t i= 0; i < dim; ++i)
            v[i]-= (Block[i]*v[dim] + Block2[i]*v[dim + 1])/DiagA[i];
    }
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
void PVankaSmootherCL::GSSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const
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
void PVankaSmootherCL::LRSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const
{
    DMatrixCL<double> M(SetupLocalProblem(id, NodeListVel, A, B));
    gauss_pivot (M, v);
}

template <typename Mat, typename Vec>
void
PVankaSmootherCL::Apply(const Mat& A, const Mat& B, const Mat& BT, const Mat&, Vec& x, Vec& y, const Vec& f, const Vec& g) const
{
    NodeListVelT NodeListVel;
    NodeListVel.reserve(200);
    size_t dim;
    const IdxDescCL* p1x= 0;
    if (idx_)
        for (MLIdxDescCL::const_iterator it= idx_->begin(); it != idx_->end(); ++it)
            if (it->IsExtended() && it->NumUnknowns() == g.size())
                p1x= &*it;

    // Cycle by pressure variables
    for (size_t id= 0; (!p1x && id < g.size()) || (p1x && id < p1x->GetXidx().GetNumUnknownsStdFE()); ++id ) {
        size_t id2( p1x ? p1x->GetXidx()[id] : NoIdx); // possible extended pressure unknown

        dim = B.row_beg(id+1) - B.row_beg(id);
        Vec v(dim + (id2 == NoIdx ? 1 : 2));

        NodeListVel.resize( dim);

        // copy column indices of nonzeros from row id of B to NodeListVel
        std::memcpy( &NodeListVel[0], B.raw_col() + B.row_beg( id), dim*sizeof(size_t));

        // v = resid of (A & B^T \\ B & 0) * (x \\ y) - (f \\ g) for the rows NodeListVel, id
        for ( size_t i=0; i<dim; ++i ) {
            const int irow = NodeListVel[i];
            v[i]= mul_row( A, x, irow) + mul_row( BT, y, irow) - f[irow];
        }
        v[dim]= mul_row( B, x, id) - g[id];
        if (id2!=NoIdx)
            v[dim +  1]= mul_row( B, x, id2) - g[id2];
            
        // smoothing
        switch (vanka_method_) {
            case 0 : {
                DiagSmoother(id, NodeListVel, A, B, v, id2);
            } break;
            case 2: {
                LRSmoother  (id, NodeListVel, A, B, v);
            } break;
            case 3: case 4 : {
                GSSmoother  (id, NodeListVel, A, B, v);
            } break;
        }

        // Cycle by local unknowns: correction of the approximation x
        for (size_t i= 0; i < dim; ++i) {
            x[NodeListVel[i]]-= tau_ * v[i];
        }
        y[id]-= tau_ * v[dim];
        if (id2 != NoIdx)
            y[id2]-= tau_ * v[dim + 1];
    }
}

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

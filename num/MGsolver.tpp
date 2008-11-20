/// \file
/// \brief classes that constitute the poisson/stokes-problem with MG-solver
/// \author Sven Gross, Joerg Grande, Volker Reichelt, Maxim Larin, Patrick Esser, IGPM

namespace DROPS {

template <class SmootherCL, class DirectSolverCL>
void
MGM(const MLMatrixCL::const_iterator& begin, const MLMatrixCL::const_iterator& fine, const MLMatrixCL::const_iterator& P,
     VectorCL& x, const VectorCL& b, const SmootherCL& Smoother, Uint smoothSteps,
     DirectSolverCL& Solver, int numLevel, int numUnknDirect)
{
    MLMatrixCL::const_iterator coarse= fine;
    MLMatrixCL::const_iterator coarseP= P;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : x.size() <= static_cast<Uint>(numUnknDirect) )
       || fine==begin)
    { // use direct solver
        Solver.Solve( *fine, x, b);
/*        std::cerr << "MGM: direct solver: iterations: " << Solver.GetIter()
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

template<class SmootherCL, class DirectSolverCL>
void MG(const MLMatrixCL& MGData, const MLMatrixCL& Prolong, const SmootherCL& smoother, DirectSolverCL& solver,
        VectorCL& x, const VectorCL& b, int& maxiter, double& tol, const bool residerr, Uint sm, int lvl)
{
    MLMatrixCL::const_iterator finest= --MGData.end();
    MLMatrixCL::const_iterator finestProlong= --Prolong.end();
    double resid= -1;
    double old_resid;
    VectorCL tmp;
    if (residerr == true) {
        resid= norm( b - *finest * x);
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
        MGM( MGData.begin(), finest, finestProlong, x, b, smoother, sm, solver, lvl, -1);
        if (residerr == true) {
            old_resid= resid;
            resid= norm( b - *finest * x);
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
void StokesMGM( const MLMatrixCL::const_iterator& beginA,  const MLMatrixCL::const_iterator& fineA,
                const MLMatrixCL::const_iterator& fineB,   const MLMatrixCL::const_iterator& fineBT, 
                const MLMatrixCL::const_iterator& fineprM, const MLMatrixCL::const_iterator& PVel,
                const MLMatrixCL::const_iterator& PPr, VectorCL& u, VectorCL& p, const VectorCL& b,
                const VectorCL& c, const StokesSmootherCL& Smoother, Uint smoothSteps, Uint cycleSteps,
                StokesDirectSolverCL& Solver, int numLevel, int numUnknDirect)
{
    MLMatrixCL::const_iterator coarseA   = fineA;
    MLMatrixCL::const_iterator coarseB   = fineB;
    MLMatrixCL::const_iterator coarseBT  = fineBT;
    MLMatrixCL::const_iterator coarseprM = fineprM;
    MLMatrixCL::const_iterator coarsePVel= PVel;
    MLMatrixCL::const_iterator coarsePPr = PPr;

    if(  ( numLevel==-1      ? false : numLevel==0 )
       ||( numUnknDirect==-1 ? false : u.size() <= static_cast<Uint>(numUnknDirect) )
       || fineA==beginA)
    {
        // use direct solver
        std::cout << "P2P1:StokesMGM: use direct solver " << std::endl;
        Solver.Solve( *fineA, *fineB, u, p, b, c);
        //std::cerr << "P2P1:StokesMGM: direct solver: iter: " << Solver.GetIter() << "\tresid: " << Solver.GetResid() << std::endl;
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

template <typename Mat, typename Vec>
void PVankaSmootherCL::DiagSmoother(Uint id, NodeListVelT& NodeListVel, const Mat& A, const Mat& B, Vec& v) const
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

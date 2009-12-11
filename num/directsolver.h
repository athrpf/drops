/// \file directsolver.h
/// \brief direct solvers
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande; SC RWTH Aachen:

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

#ifndef DROPS_DIRECTSOLVER_H
#define DROPS_DIRECTSOLVER_H

#include "num/spmat.h"
#include "cholmod.h"
#include "umfpack.h"

namespace DROPS {

/*******************************************************************
*   D I R E C T N O N S Y M M S O L V E R   C L                    *
*******************************************************************/
/// \brief Use direct solvers for Ax=b with sparse nonsymmetric matrix A.
/** This class offers an interface to the UMFPACK package from
    http://www.cise.ufl.edu/research/sparse/                       */
/*******************************************************************
*   D I R E C T N O N S Y M M S O L V E R   C L                    *
*******************************************************************/
class DirectNonSymmSolverCL {

private:
    double Info_ [UMFPACK_INFO],        ///< stores status-infos
           Control_ [UMFPACK_CONTROL];  ///< control parameters for UMFPACK
    void   *Symbolic_,                  ///< stores symbolic factorization
           *Numeric_;                   ///< stores numeric factorization
    int    *Apt,                        ///< matrix A row indices
           *Ait;                        ///< matrix A column pointer
    double *Axt;                        ///< matrix A nonzeros
    Uint   num_rows_,                   ///< # rows
           num_cols_;                   ///< # columns

public:

/// store and factorize the matrix A
DirectNonSymmSolverCL(const MatrixCL& A)
: Apt(0), Ait(0), Axt(0)
{
    // get the default control parameters
    umfpack_di_defaults(Control_);

    // change the default print level, uncomment the following line for verbose output
    //Control_ [UMFPACK_PRL] = 6;

    // print the license agreement
    umfpack_di_report_status (Control_, UMFPACK_OK);

    // print the control parameters
    umfpack_di_report_control (Control_);

    Update(A);
}

~DirectNonSymmSolverCL()
{
    delete [] Apt;
    delete [] Ait;
    delete [] Axt;
    umfpack_di_free_symbolic (&Symbolic_);
    umfpack_di_free_numeric (&Numeric_);
}
/// delete stored matrix and store/factorize the matrix A
void Update(const MatrixCL& A)
{
    int status;

    // DROPS: CCS, UMFPACK CRS
    int* Ap    = new int[A.num_rows()+1];
    int* Ai    = new int[A.num_nonzeros()];
    double* Ax = new double[A.num_nonzeros()];

    std::copy(A.raw_col(), A.raw_col()+A.num_nonzeros(), Ai);
    std::copy(A.raw_row(), A.raw_row()+A.num_rows()+1  , Ap);
    std::copy(A.raw_val(), A.raw_val()+A.num_nonzeros(), Ax);

    delete [] Apt;
    delete [] Ait;
    delete [] Axt;

    Apt = new int    [A.num_rows()+1];
    Ait = new int    [A.num_nonzeros()];
    Axt = new double [A.num_nonzeros()];

    status = umfpack_di_transpose(A.num_cols(), A.num_rows(), Ap, Ai, Ax, 0, 0, Apt, Ait, Axt);
    if (status < 0)
    {
        umfpack_di_report_info (Control_, Info_);
        umfpack_di_report_status (Control_, status);
        throw DROPSErrCL("umfpack_di_transpose failed") ;
    }

    delete [] Ap;
    delete [] Ai;
    delete [] Ax;

    umfpack_di_report_matrix (A.num_rows(), A.num_cols(), Apt, Ait, Axt, 1, Control_);

    // symbolic factorization
    status = umfpack_di_symbolic (A.num_rows(), A.num_cols(), Apt, Ait, Axt, &Symbolic_, Control_, Info_) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control_, Info_);
        umfpack_di_report_status (Control_, status);
        throw DROPSErrCL("umfpack_di_symbolic failed");
    }

    //print the symbolic factorization
    umfpack_di_report_symbolic (Symbolic_, Control_);

    // numeric factorization
    status = umfpack_di_numeric (Apt, Ait, Axt, Symbolic_, &Numeric_, Control_, Info_);
    if (status < 0)
    {
        umfpack_di_report_info (Control_, Info_);
        umfpack_di_report_status (Control_, status);
        throw DROPSErrCL("umfpack_di_numeric failed") ;
    }

    // print the numeric factorization
    umfpack_di_report_numeric (Numeric_, Control_);

    num_cols_ = A.num_cols();
    num_rows_ = A.num_rows();
}

/// solve Ax=b
void Solve(const MatrixCL, VectorCL &x, const VectorCL& b)
{
    int status;
    if (b.size() != num_rows_)
        throw DROPSErrCL("DirectNonSymmSolverCL::Solve: incompatible dimensions");

    status = umfpack_di_solve (UMFPACK_A, Apt, Ait, Axt, Addr(x), Addr(b), Numeric_, Control_, Info_);
    umfpack_di_report_info (Control_, Info_);
    umfpack_di_report_status (Control_, status);
    if (status < 0)
    {
        throw DROPSErrCL("umfpack_di_solve failed");
    }
    umfpack_di_report_vector (num_cols_, Addr(x), Control_);
}
};

/*******************************************************************
*   D I R E C T S Y M M S O L V E R   C L                          *
*******************************************************************/
/// \brief Use direct solvers for Ax=b with sparse symmetric matrix A.
/** This class offers an interface to the CHOLMOD package from
    http://www.cise.ufl.edu/research/sparse/                       */
/*******************************************************************
*   D I R E C T       S Y M M S O L V E R   C L                    *
*******************************************************************/
class DirectSymmSolverCL {

private:
    cholmod_sparse *A_;                 ///< sparse matrix A
    cholmod_factor *L_;                 ///< factor of A
    cholmod_common c_;                  ///< CHOLMOD parameter

public:

/// store and factorize the matrix A
DirectSymmSolverCL(const MatrixCL& A) :
A_(0), L_(0)
{
    cholmod_l_start(&c_);
    Update(A);
}

/// delete stored matrix and store/factorize the matrix A
void Update(const MatrixCL& A)
{
    cholmod_l_free_sparse (&A_, &c_);
    A_= cholmod_l_allocate_sparse (A.num_rows(), A.num_cols(), A.num_nonzeros(), /*sorted*/ true,
           /*packed*/ true, /*upper left block used*/ 1, /*pattern*/ CHOLMOD_REAL, &c_);

    // MatrixCL -> cholmod_sparse
    std::copy(A.raw_col(), A.raw_col()+A.num_nonzeros(), static_cast<size_t *>(A_->i));
    std::copy(A.raw_row(), A.raw_row()+A.num_rows()+1  , static_cast<size_t *>(A_->p));
    std::copy(A.raw_val(), A.raw_val()+A.num_nonzeros(), static_cast<double *>(A_->x));

   //check A
    cholmod_l_check_sparse(A_, &c_);
//    char dummy[]= "A";
//    cholmod_print_sparse (A_, dummy, &c_);

    L_ = cholmod_l_analyze (A_, &c_);
    cholmod_l_factorize (A_, L_, &c_);
}

/// solve Ax=b
void Solve(const MatrixCL, VectorCL &x, const VectorCL& b)
{
    cholmod_dense  *b_;
    cholmod_dense  *x_;
    b_= cholmod_l_allocate_dense (b.size(), 1, /*leading dim.*/ b.size(), /*pattern*/ CHOLMOD_REAL, &c_);

    // VectorCL -> cholmod_dense
    std::copy(Addr(b), Addr(b)+b.size(), static_cast<double *>(b_->x));

    // check b
    cholmod_check_dense(b_, &c_);
//    char dummy[]= "b";
//    cholmod_print_dense (b_, dummy, &c_);

    // solve Ax=b
    x_ = cholmod_l_solve (CHOLMOD_A, L_, b_, &c_);

    // cholmod_dense -> VectorCL
    x = VectorCL(static_cast<double *>(x_->x), b.size());
    cholmod_l_free_dense  (&b_, &c_);
    cholmod_l_free_dense  (&x_, &c_);
}

~DirectSymmSolverCL()
{
    cholmod_l_free_factor (&L_, &c_);
    cholmod_l_free_sparse (&A_, &c_);
    cholmod_l_finish(&c_);
}

};

}//end of namespace DROPS

#endif

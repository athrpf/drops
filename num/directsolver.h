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
    double  Info_ [UMFPACK_INFO],        ///< stores status-infos
            Control_ [UMFPACK_CONTROL];  ///< control parameters for UMFPACK
    void    *Symbolic_,                  ///< stores symbolic factorization
            *Numeric_;                   ///< stores numeric factorization
    UF_long *Apt,                        ///< matrix A row indices
            *Ait;                        ///< matrix A column pointer
    double  *Axt;                        ///< matrix A nonzeros
    size_t  num_rows_,                   ///< # rows
            num_cols_;                   ///< # columns

public:

/// store and factorize the matrix A
DirectNonSymmSolverCL(const MatrixCL& A)
: Apt(0), Ait(0), Axt(0)
{
    // get the default control parameters
    umfpack_dl_defaults(Control_);

    // change the default print level, uncomment the following line for verbose output
    //Control_ [UMFPACK_PRL] = 6;

    // print the license agreement
    umfpack_dl_report_status (Control_, UMFPACK_OK);

    // print the control parameters
    umfpack_dl_report_control (Control_);

    Update(A);
}

~DirectNonSymmSolverCL()
{
    delete [] Apt;
    delete [] Ait;
    delete [] Axt;
    umfpack_dl_free_symbolic (&Symbolic_);
    umfpack_dl_free_numeric (&Numeric_);
}
/// delete stored matrix and store/factorize the matrix A
void Update(const MatrixCL& A)
{
    int status;

    delete [] Apt;
    delete [] Ait;
    delete [] Axt;

    Apt = new UF_long [A.num_rows()+1];
    Ait = new UF_long [A.num_nonzeros()];
    Axt = new double  [A.num_nonzeros()];

    std::copy(A.raw_col(), A.raw_col()+A.num_nonzeros(), Ait);
    std::copy(A.raw_row(), A.raw_row()+A.num_rows()+1  , Apt);
    std::copy(A.raw_val(), A.raw_val()+A.num_nonzeros(), Axt);

    //print the matrix
    umfpack_dl_report_matrix (A.num_rows(), A.num_cols(), Apt, Ait, Axt, 1, Control_);

    // symbolic factorization
    status = umfpack_dl_symbolic (A.num_rows(), A.num_cols(), Apt, Ait, Axt, &Symbolic_, Control_, Info_) ;
    if (status < 0)
    {
        Control_ [UMFPACK_PRL] = 6;
        umfpack_dl_report_info   (Control_, Info_);
        umfpack_dl_report_status (Control_, status);
        throw DROPSErrCL("umfpack_dl_symbolic failed");
    }

    //print the symbolic factorization
    umfpack_dl_report_symbolic (Symbolic_, Control_);

    // numeric factorization
    status = umfpack_dl_numeric (Apt, Ait, Axt, Symbolic_, &Numeric_, Control_, Info_);
    if (status < 0)
    {
        Control_ [UMFPACK_PRL] = 6;
        umfpack_dl_report_info   (Control_, Info_);
        umfpack_dl_report_status (Control_, status);
        throw DROPSErrCL("umfpack_dl_numeric failed") ;
    }

    // print the numeric factorization
    umfpack_dl_report_numeric (Numeric_, Control_);

    num_cols_ = A.num_cols();
    num_rows_ = A.num_rows();
}

/// solve Ax=b, NOTE: DROPS::CRS, UMFPACK: CCS
void Solve(const MatrixCL, VectorCL &x, const VectorCL& b)
{
    int status;
    if (b.size() != num_cols_)
        throw DROPSErrCL("DirectNonSymmSolverCL::Solve: incompatible dimensions");

    status = umfpack_dl_solve (UMFPACK_At, Apt, Ait, Axt, Addr(x), Addr(b), Numeric_, Control_, Info_);
    if (status < 0)
    {
        Control_ [UMFPACK_PRL] = 6;
        umfpack_dl_report_info   (Control_, Info_);
        umfpack_dl_report_status (Control_, status);
        throw DROPSErrCL("umfpack_dl_solve failed");
    }
    umfpack_dl_report_vector (num_rows_, Addr(x), Control_);
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

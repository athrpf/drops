/// \file hypre.cpp
/// \brief Interface to HYPRE (high performance preconditioners), Lawrence Livermore Nat. Lab. 
/// \author LNM RWTH Aachen: Sven Gross; SC RWTH Aachen: Oliver Fortmeier

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

#ifdef _HYPRE
#include "num/hypre.h"

namespace DROPS {

void HypreIndexCL::Clear()
{
    idx_beg_.clear();
    global_idx_.resize(0);
    exclusive_idx_.clear();
    local_idx_.clear();
    exclusive_global_idx_.resize(0);
}

/** To generate a HYPRE index, a global indexing (int) is performed. This involves a
    collective communication and a neighbor communication.*/
void HypreIndexCL::Generate()
{
    const size_t numUnk= idx_.NumUnknowns();
    num_exclusive_= numUnk;
    const ExchangeCL& ex= idx_.GetEx();
    num_exclusive_= ex.GetNumExclusive();

    Clear();
    idx_beg_.resize( ProcCL::Size()+1);
    global_idx_.resize( numUnk, -1);
    exclusive_idx_.resize( numUnk, -1);
    local_idx_.resize(  numUnk, -1);
    exclusive_global_idx_.resize( num_exclusive_);
    
    // assign each exclusive index a consecutive number
    int gidx=0;
    for ( size_t i=0; i<numUnk; ++i){
        if ( ex.IsExclusive( i)){
            global_idx_[i]= gidx++;
            if (gidx<0) throw DROPSErrCL("HypreIndexCL::Generate: int is to small");
        }
    }
    if (gidx!=(int)num_exclusive_){ 
        throw DROPSErrCL("HypreIndexCL::Generate: Number of exclusive dof does not match");
    }

    // Distribute information among processes
    std::vector<int> idx_rcv(ProcCL::Size());
    ProcCL::Gather( gidx, Addr(idx_rcv));

    // Fill beginning indices
    idx_beg_[0]=0;
    for ( int proc=0; proc<ProcCL::Size(); ++proc){
        idx_beg_[proc+1]= idx_beg_[proc]+idx_rcv[proc];
    }

    // Generate a global numbering
    for ( size_t i=0; i<numUnk; ++i){
        if (global_idx_[i]==-1)
            global_idx_[i]=0;
        else
            global_idx_[i] += GetLower();
    }
    ex.Accumulate( global_idx_);

    // Generate mapping global index to local index and local exclusive index
    gidx=0;
    for ( size_t i=0; i<numUnk; ++i){
        if ( IsLocExclusive(i)){
            local_idx_[    global_idx_[i]-GetLower()]= (int) i;
            exclusive_idx_[global_idx_[i]-GetLower()]= gidx++;
        }
    }

#if DROPSDebugC & DebugNumericC
    // Check if mapping is correct
    for ( size_t i=0; i<numUnk; ++i){
        if ( IsLocExclusive(i)){            // i is exclusive
            if ( GetLocalIdx( GetGlobalIdx(i)) != (int) i)
                throw DROPSErrCL("HypreIndexCL::Generate: Mapping global <-> local does not match for an exclusive index");
            if ( GetExclusiveNum(GetGlobalIdx(i))<0 ||  GetExclusiveNum(GetGlobalIdx(i))>=(int)num_exclusive_)
                throw DROPSErrCL("HypreIndexCL::Generate: Mapping global <-> local exclusive does not match for an exclusive index");
        }
        else{                               // i is located on another process
            if ( GetLocalIdx( GetGlobalIdx(i)) != -1)
                throw DROPSErrCL("HypreIndexCL::Generate: Mapping global <-> local does not match for a remote index");
            if ( GetExclusiveNum(GetGlobalIdx(i)) != -1)
                throw DROPSErrCL("HypreIndexCL::Generate: Mapping global <-> local exclusive does not match for a remote index");
        }
    }
#endif

    // store global numbering of exclusive dof
    gidx=0;
    for ( int i=0; i<(int)global_idx_.size(); ++i)
        if ( IsGlobExclusive( global_idx_[i]))
            exclusive_global_idx_[gidx++]= global_idx_[i];
}

void HypreMatrixCL::CreateSendMPIType()
{
    if (sendMPItype_ != ProcCL::NullDataType)
        return;

    NonZeroST tmp;
    const int count               = 3;
    const int block[3]            = {1,1,1};
    const ProcCL::AintT sAddr     = ProcCL::Get_address(&tmp);
    const ProcCL::DatatypeT typ[4]= { ProcCL::MPI_TT<int>::dtype, ProcCL::MPI_TT<int>::dtype,
                                      ProcCL::MPI_TT<double>::dtype};
    ProcCL::AintT displace[3];

    displace[0]= ProcCL::Get_address(&tmp.row)  - sAddr;
    displace[1]= ProcCL::Get_address(&tmp.col)  - sAddr;
    displace[2]= ProcCL::Get_address(&tmp.value)- sAddr;

    sendMPItype_ = ProcCL::CreateStruct(count, block, displace, typ);
    ProcCL::Commit(sendMPItype_);
}

void HypreMatrixCL::Clear()
{
    ncols_.resize(0);
    nOffCol_.resize( 0);
    nMyCol_.resize( 0);
    cols_.resize(0);
    vals_.resize(0);
}

/** Calling of HYPRE functions to set up the HYPRE data structures and the MPI type
    for sending a non-zero.*/
void HypreMatrixCL::Init()
{
    HYPRE_IJMatrixCreate( ProcCL::GetComm(), rowidx_.GetLower(), rowidx_.GetUpper(), colidx_.GetLower(), colidx_.GetUpper(), &ijMat_);
    HYPRE_IJMatrixSetObjectType( ijMat_, HYPRE_PARCSR);
    CreateSendMPIType();
}


/** This function distributes the given matrix row-wise among the processes. To this
    end, the HYPRE index object corresoponding to the rows is used. This function
    takes care of sending the non-zero entries among the processes and calling the
    HYPRE functions to generate the matrix in the HYPRE format.*/
void HypreMatrixCL::Generate()
{
    Clear();
    typedef std::map< int, std::vector<NonZeroST> > SendBufT;
    SendBufT sendBuf;
    // send at least each neighbor a message
    for ( ExchangeCL::ProcNumCT::const_iterator neigh= rowidx_.GetEx().GetNeighbors().begin(); neigh!=rowidx_.GetEx().GetNeighbors().end(); ++neigh){
        sendBuf[*neigh].push_back( NonZeroST());
    }
    // iterate over matrix and fill sendbuffer if row index is not exclusive
    for ( size_t i=0; i<M_.num_rows(); ++i){
        if ( !rowidx_.IsLocExclusive( i)){
            int exProc= rowidx_.GetLocExclusiveProc(i);
            for ( size_t nz=M_.row_beg( i); nz<M_.row_beg( i+1); ++nz){
                sendBuf[exProc].push_back( NonZeroST( rowidx_.GetGlobalIdx( i), colidx_.GetGlobalIdx( M_.col_ind(nz)), M_.val( nz)));
            }
        }
    }

    // initiate sending
    std::vector< ProcCL::RequestT> sendReq;
    for ( SendBufT::const_iterator it=sendBuf.begin(); it!=sendBuf.end(); ++it){
        sendReq.push_back( ProcCL::Isend( Addr(it->second), it->second.size(), sendMPItype_, it->first, tag_));
    }

    typedef std::vector< std::vector<int> >    ColT;        // Store column indices of each exclusive row
    ColT colind( rowidx_.GetNumExclusive());
    typedef std::vector< std::vector<double> > ValT;        // Store corresponding non-zeros
    ValT values( rowidx_.GetNumExclusive());

    // fill colind and values with local elements
    size_t localrowpos=0;
    for ( size_t i=0; i<M_.num_rows(); ++i){
        if ( rowidx_.IsLocExclusive( i)){
            for ( size_t nz=M_.row_beg( i); nz<M_.row_beg( i+1); ++nz){
                colind[localrowpos].push_back( colidx_.GetGlobalIdx(M_.col_ind( nz)));
                values[localrowpos].push_back( M_.val( nz));
            }
            ++localrowpos;
        }
    }

    // receive elements
    std::vector<NonZeroST> recvBuf;
    for ( ExchangeCL::ProcNumCT::const_iterator neigh= rowidx_.GetEx().GetNeighbors().begin(); neigh!=rowidx_.GetEx().GetNeighbors().end(); ++neigh){
        recvBuf.resize( ProcCL::GetMessageLength( *neigh, tag_, sendMPItype_));
        ProcCL::Recv( Addr(recvBuf), recvBuf.size(), sendMPItype_, *neigh, tag_);
        for ( size_t i=0; i<recvBuf.size(); ++i){
            if ( !recvBuf[i].IsDummy()){
                const int localRow= rowidx_.GetExclusiveNum( recvBuf[i].row);
                Assert( localRow!=-1, DROPSErrCL("HypreMatrixCL::Generate: Received elements, I am not responsible for"), DebugParallelNumC);
                // Check if an entry for this column already exists
                std::vector<int>::iterator pos= std::find( colind[ localRow].begin(), colind[ localRow].end(), recvBuf[i].col);
                if ( pos==colind[ localRow].end()){
                    colind[ localRow].push_back( recvBuf[i].col);
                    values[ localRow].push_back( recvBuf[i].value);
                }
                else{
                    int nzPos= std::distance(colind[ localRow].begin(), pos);
                    Assert(nzPos>=0 && nzPos<(int)values[localRow].size(), DROPSErrCL("HypreMatrixCL::Generate: Non-zero does not exists"), DebugParallelNumC);
                    values[ localRow][nzPos] += recvBuf[i].value;
                }
            }
        }
    }

    // generate data format to pass to HYPRE
    ncols_.resize( rowidx_.GetNumExclusive(), 0);
    nOffCol_.resize( rowidx_.GetNumExclusive(), 0);
    nMyCol_.resize( rowidx_.GetNumExclusive(), 0);
    size_t nnz=0;
    for ( int row=0; row<rowidx_.GetNumExclusive(); ++row){
        ncols_[row] = colind[row].size();              // number of column indices per row
        nnz        += colind[row].size();              // number of non-zeros
        Assert(colind[row].size()==values[row].size(), DROPSErrCL("HypreMatrixCL::Generate: Size of values and column indices does not match"), DebugParallelNumC);
    }
    cols_.resize( nnz, 0);
    vals_.resize( nnz, 0.);
    size_t nzpos=0;
    for ( int row=0; row<rowidx_.GetNumExclusive(); ++row){
        for ( size_t col=0; col<colind[row].size(); ++col){
            const int columnindex= colind[row][col];
            vals_[nzpos  ]= values[row][col];
            cols_[nzpos++]= columnindex;
            if ( colidx_.IsGlobExclusive(columnindex))
                nMyCol_[row]++;
            else
                nOffCol_[row]++;
        }
    }

    // Call HYPRE functions to generate the parallel CRS matrix
    HYPRE_IJMatrixSetRowSizes( ijMat_, Addr(ncols_));
    HYPRE_IJMatrixSetDiagOffdSizes( ijMat_, Addr(nMyCol_), Addr(nOffCol_));
//    HYPRE_IJMatrixSetMaxOffProcElmts( ijMat_, numOffElem);

    HYPRE_IJMatrixInitialize( ijMat_);
    HYPRE_IJMatrixSetValues( ijMat_, rowidx_.GetNumExclusive(), Addr(ncols_), rowidx_.GetGlobalExList(), Addr(cols_), Addr(vals_));
    HYPRE_IJMatrixAssemble( ijMat_);
    HYPRE_IJMatrixGetObject( ijMat_, (void**) &parMat_);

    // Wait until free memory of send buffers
    ProcCL::WaitAll(sendReq);
}

void HypreVectorCL::Init()
{
    HYPRE_IJVectorCreate( ProcCL::GetComm(), idx_.GetLower(), idx_.GetUpper(), &ijVec_);
    HYPRE_IJVectorSetObjectType( ijVec_, HYPRE_PARCSR);
}

/** Pass vector entries to HYPRE. Note that this function assums the vector
    to be in accumulated format.*/
void HypreVectorCL::Generate()
{
//    HYPRE_IJVectorSetMaxOffProcElmts( ijVec_, idx_.GetOffRange() );
    VectorCL vtmp( idx_.GetNumExclusive());
    size_t pos=0;
    for ( size_t i=0; i<v_.size(); ++i)
        if ( idx_.IsLocExclusive(i))
            vtmp[pos++]= v_[i];
    HYPRE_IJVectorInitialize( ijVec_);
    HYPRE_IJVectorSetValues( ijVec_, vtmp.size(), idx_.GetGlobalExList(), Addr(vtmp));
    HYPRE_IJVectorAssemble( ijVec_);
    HYPRE_IJVectorGetObject( ijVec_, (void**) &parVec_);
}

/** Copy data from the HYPRE vector into the DROPS vector and perform a neighbor
    communication to gain an accumulated result.
*/
void HypreVectorCL::Retrieve()
{
    const VectorBaseCL<int>& global_idx= idx_.GetGlobalExListVector();
    v_= 0.;
    VectorCL vtmp( idx_.GetNumExclusive());
    HYPRE_IJVectorGetValues( ijVec_, idx_.GetNumExclusive(), idx_.GetGlobalExList(), Addr(vtmp));
    size_t pos=0;
    for ( size_t i=0; i<global_idx.size(); ++i){
        if ( idx_.IsGlobExclusive( global_idx[i])){
            v_[ idx_.GetLocalIdx( global_idx[i])]= vtmp[ pos++];
        }
    }
    idx_.GetEx().Accumulate( v_);
}

HypreVectorCL::~HypreVectorCL()
{
    HYPRE_IJVectorDestroy( ijVec_);
}

void HypreAMGSolverCL::Init()
{
    HYPRE_BoomerAMGCreate( &solver_);
}

HypreAMGSolverCL::~HypreAMGSolverCL()
{
    HYPRE_BoomerAMGDestroy( solver_);
}

void HypreAMGSolverCL::SetupAndSolve( const MatrixCL& A, VectorCL& x, const VectorCL& b, const IdxDescCL& idx) const 
{
    HypreMatrixCL hA(A, idx, idx);
    HypreVectorCL hx(x, hA.GetRowIdx());
    VectorCL acc_b= hA.GetColIdx().GetEx().GetAccumulate( b);
    HypreVectorCL hb(const_cast<VectorCL&>(acc_b), hA.GetColIdx());
    Setup( hA, hx, hb);
    Solve( hA, hx, hb);
    hx.Retrieve(); // copy solution to x
}

void HypreAMGSolverCL::Setup( const HypreMatrixCL& A, const HypreVectorCL& x, const HypreVectorCL& b) const 
{ 
    // set some solver parameters...
    HYPRE_BoomerAMGSetup( solver_, A(), b(), x()); 
}

void HypreAMGSolverCL::Solve( const HypreMatrixCL& A, HypreVectorCL& x, const HypreVectorCL& b) const 
{ 
    // set stopping criterion
    HYPRE_BoomerAMGSetTol( solver_, _tol); 
    HYPRE_BoomerAMGSetMaxIter( solver_, _maxiter);
    // solve 
    HYPRE_BoomerAMGSolve( solver_, A(), b(), x());
    // get iteration info
    HYPRE_BoomerAMGGetNumIterations( solver_, &_iter); 
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm( solver_, &_res); 
}

}    // end of namespace DROPS
#endif  // _HYPRE

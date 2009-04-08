#include "num/spmat.h"
#include "misc/problem.h"
#include <iostream>
#include <fstream>

using namespace DROPS;

void
print_frobeniusnorm (const MatrixCL& M)
{
    double F= 0.0, FL= 0.0, FU= 0.0, FD= 0.0;
    for (size_t i= 0; i< M.num_rows(); ++i)
        for (size_t nz= M.row_beg( i); nz < M.row_beg( i + 1); ++nz) {
            const double vsq( std::pow( M.val( nz), 2));
            F+= vsq;
            if ( M.col_ind( nz) < i) FL+= vsq;
            if ( M.col_ind( nz) > i) FU+= vsq;
            if ( M.col_ind( nz) == i) FD+= vsq;
        }
    std::cout << "frobeniusnorm: F: " << std::sqrt( F) << '\t'
        << "FL: " << std::sqrt( FL) << '\t'
        << "FU: " << std::sqrt( FU) << '\t'
        << "FD: " << std::sqrt( FD) << '\n';
}

// Computes B with B_ij= A_ij - max( A_ij, A_ji).
void
non_symmetric_part (const MatrixCL& M, MatrixCL& NP)
{
    SparseMatBuilderCL<double> N( &NP, M.num_rows(), M.num_cols());
    for (size_t i= 0; i< M.num_rows(); ++i) {
        for (size_t nz= M.row_beg( i); nz < M.row_beg( i + 1); ++nz) {
            const double Mij= M.val( nz);
            const double Mji= M( M.col_ind( nz), i);
            if ( Mij < Mji ) N( i, M.col_ind( nz))= Mij - Mji;
        }
    }
    N.Build();
    std::cout << "non_symmetric_part: nonzeros: " << NP.num_nonzeros () << '\n';
}

void
sort_row_entries (MatrixCL& M)
{
    // Sort the indices in each row by increasing value of the entries.
    // Note, that this is dangerous as most parts of Drops assume ordering by column-index.
    // It is required for the renumbering algo.
    typedef std::pair<double, size_t> PT;
    for (size_t r= 0; r < M.num_rows(); ++r) {
        std::vector<PT> pv( M.row_beg( r + 1) - M.row_beg( r));
        for (size_t i= M.row_beg( r), j= 0; i < M.row_beg( r + 1); ++i, ++j)
            pv[j]= std::make_pair( M.val( i), M.col_ind( i));
        std::sort( pv.begin(), pv.end(), less1st<PT>());
        for (size_t i= M.row_beg( r), j= 0; i < M.row_beg( r + 1); ++i, ++j) {
            M.raw_col()[i]= pv[j].second;
            M.raw_val()[i]= pv[j].first;
        }        
    }
}

VectorCL
weight_in (const MatrixCL& M)
{
    VectorCL w( 1.0, M.num_rows());
    return M*w;
}

VectorCL
weight_out (const MatrixCL& M)
{
    VectorCL w( 1.0, M.num_rows());
    return transp_mul( M, w);
}

PermutationT
sort_vertices_by_weight (const VectorCL& w)
{
    typedef std::pair<double, size_t> PT;

    std::vector<PT> pv( w.size());
    for (size_t r= 0; r < w.size(); ++r)
        pv[r]= std::make_pair( w[r], r);
    std::sort( pv.begin(), pv.end(), less1st<PT>());

    PermutationT p( w.size());
    std::transform( pv.begin(), pv.end(), p.begin(), select2nd<PT>());
//    for (size_t r= 0; r < w.size(); ++r)
//            p[r]= pv[r].second;
    return p;
}


void
depth_first (IdxT v, const MatrixCL& M, PermutationT& p, IdxT& idx, size_t& c)
{
    IdxT w;
    --p[v];
    for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
        if (M.val( e) == 0.0 || (w= M.col_ind( e)) == v)
            continue;
        if (p[w] == NoIdx)
            depth_first( w, M, p, idx, c);
        if (p[w] == NoIdx - 1) {
            ++c;
//            std::cout << "Cycle detected at vertex " << w
//                << ". Edge (" << v << ", " << w << ") ignored.\n";
        }
    }
    p[v]= idx++;
}

// M is interpreted as digraph. Vertices: x_i, Edges: (i,j) iff M_ji != 0
// Diagonal entries are ignored (Cycles of length 1).
// If M is acyclic, p contains a permutation of 0..M.num_rows(), such
// that i<j for all edges (i,j); that means, M is lower triangular.
// If M contains cycles, these are broken by ignoring the closing edge
// in the depth-first-search.
//
// NoIdx is used as "not visited"-mark in p.
// NoIdx - 1 is used as "visited, but not finnished"-mark in p.
// If p.size() >= NoIdx - 1 an error is raised, but probably DROPS would
// segfault before that anyways (The corresponding matrix would be huge).
void
topo_sort (const MatrixCL& M, PermutationT& p)
{
    if (M.num_rows() >= NoIdx - 1)
        throw DROPSErrCL( "topo_sort: Graph is too large.\n");

    VectorCL w( weight_in( M) - weight_out( M));
std::ofstream ow( "weights.dat");
ow << w;
    PermutationT pw( sort_vertices_by_weight ( w));

    size_t c= 0;
    IdxT idx= 0;
    p.assign( M.num_rows(), NoIdx);
    for (IdxT v= 0; v < M.num_rows(); ++v)
        if (p[pw[v]] == NoIdx) depth_first( pw[v], M, p, idx, c);

    if (idx != M.num_rows())
        throw DROPSErrCL( "topo_sort: Sort failed.\n");

    std::cout << "topo_sort: " << c << " cycles detected.\n";
}


void
bey_depth_first (IdxT v, const MatrixCL& M, const MatrixCL& Mt, PermutationT& p, IdxT& idx)
{
    IdxT w;

    // Check, that all predeccessors are numbered.
    for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
        if (M.val( e) == 0.0 || M.col_ind( e) == v)
            continue;
        if (p[M.col_ind( e)] == NoIdx) return;
    }

    p[v]= idx++;

    // Try to number the successors.
    for (size_t e= Mt.row_beg( v); e < Mt.row_beg( v + 1); ++e) {
        if (Mt.val( e) == 0.0 || (w= Mt.col_ind( e)) == v)
            continue;
        if (p[w] == NoIdx)
            bey_depth_first( w, M, Mt, p, idx);
    }
}

void
bey (const MatrixCL& M, PermutationT& p)
{
    if (M.num_rows() >= NoIdx - 1)
        throw DROPSErrCL( "bey: Graph is too large.\n");

    MatrixCL Mt;
    transpose( M, Mt);
    IdxT idx= 0;
    p.assign( M.num_rows(), NoIdx);
    for (IdxT v= 0; v < M.num_rows(); ++v)
        if (p[v] == NoIdx) bey_depth_first( v, M, Mt, p, idx);

    std::cout << "bey: Could number " << idx << " dofs topologically.\n";
    for (IdxT v= 0; v < M.num_rows(); ++v)
        if (p[v] == NoIdx) p[v]= idx++;

    if (idx != M.num_rows())
        throw DROPSErrCL( "topo_sort: Sort failed.\n");
}


void TestPermutation()
{
    MatrixCL M;
    MatrixBuilderCL Mb( &M, 3, 3);
    Mb( 0, 0)= 0.0;
    Mb( 1, 1)= 0.0;
    Mb( 2, 2)= 0.0;

    Mb( 0, 1)= 1.0;
    Mb( 0, 2)= 1.0;

    Mb( 1, 0)= 0.0;
    Mb( 1, 2)= 1.0;
    Mb.Build();
    std::cout << "M (the digraph):\n" << M << '\n';

    PermutationT p;
    topo_sort ( M, p);
    seq_out( p.begin(), p.end(), std::cout);

    VectorCL v( 3);
    v[0]= 0.; v[1]= 1.; v[2]= 2.;
    permute_Vector( v, p);
    std::cout << v << '\n';

    MatrixCL m;
    MatrixBuilderCL mb( &m, 3, 3);
    mb( 0, 0)= 1.0;
    mb( 0, 1)= 2.0;
    mb( 0, 2)= 3.0;
    mb( 1, 0)= 4.0;
    mb( 1, 1)= 5.0;
    mb( 1, 2)= 6.0;
    mb( 2, 0)= 7.0;
    mb( 2, 1)= 8.0;
    mb( 2, 2)= 9.0;
    mb.Build();
    MatrixCL m2( m), m3( m);
    
    std::cout << "Matrix:\n" << m << '\n' << "Permutation:\n";
    seq_out( p.begin(), p.end(), std::cout);

    m.permute_columns( p);
    std::cout << "Columns permuted:\n" << m << '\n';

    m2.permute_rows( p);
    std::cout << "Rows permuted:\n" << m2 << '\n';

    m3.permute_columns( p);
    m3.permute_rows( p);
    std::cout << "Columns and rows permuted:\n" << m3 << '\n';

    
    print_frobeniusnorm( M);
    M.permute_rows( p);
    std::cout << "M, rows permuted:\n" << M << '\n';
    M.permute_columns( p);
    std::cout << "M permuted:\n" << M << '\n';
    print_frobeniusnorm( M);
}

void
depth_first_strategy()
{
    std::ifstream in( "normal.dat");
    if (!in) throw DROPSErrCL( "Datei nicht vorhanden.\n");
    MatrixCL M;
    in >> M;
    // std::ofstream of( "normal2out.dat");
    // of << M;
    MatrixCL N;
    non_symmetric_part( M, N);
    sort_row_entries ( N);
    std::ofstream ofNP( "normalNP.dat");
    ofNP << N;
    print_frobeniusnorm( M);
    PermutationT p;
    topo_sort ( N, p);
    // seq_out( p.begin(), p.end(), std::cout);
    M.permute_rows( p);
    M.permute_columns( p);
    print_frobeniusnorm( M);
    std::ofstream of( "normalNPperm.dat");
    of << M;
}

void
bey_strategy()
{
    std::ifstream in( "normal.dat");
    if (!in) throw DROPSErrCL( "Datei nicht vorhanden.\n");
    MatrixCL M;
    in >> M;
    // std::ofstream of( "normal2out.dat");
    // of << M;
    MatrixCL N;
    non_symmetric_part( M, N);
    sort_row_entries ( N);
    std::ofstream ofNP( "normalNP.dat");
    ofNP << N;
    print_frobeniusnorm( N);
    PermutationT p;
    bey( N, p);
    // seq_out( p.begin(), p.end(), std::cout);
    M.permute_rows( p);
    M.permute_columns( p);
    print_frobeniusnorm( M);
    std::ofstream of( "normalNPperm.dat");
    of << M;
}

void
Test_rcm()
{
    MatrixCL M;
    MatrixBuilderCL Mb( &M, 3, 3);
    Mb( 0, 0)= 0.0;
    Mb( 1, 1)= 0.0;
    Mb( 2, 2)= 0.0;

    Mb( 0, 1)= 1.0;
    Mb( 0, 2)= 1.0;

    Mb( 1, 0)= 0.0;
    Mb( 1, 2)= 1.0;
    Mb.Build();
    std::cout << "M (the digraph):\n" << M << '\n';
    print_frobeniusnorm( M);

    PermutationT p;
    reverse_cuthill_mckee( M, p);
    M.permute_rows( p);
    M.permute_columns( p);
    std::cout << "M permuted:\n" << M << '\n';
    print_frobeniusnorm( M);
}

int main ()
{
try {
//    TestPermutation();
//    depth_first_strategy();
//    bey_strategy();
//     std::ifstream in( "normal.dat");
//     if (!in) throw DROPSErrCL( "Datei nicht vorhanden.\n");
//     MatrixCL M;
//     in >> M;
//     print_frobeniusnorm( M);
//     PermutationT p;
//     reverse_cuthill_mckee( M, p);
//     M.permute_rows( p);
//     M.permute_columns( p);
//     print_frobeniusnorm( M);
//     std::ofstream of( "normalperm.dat");
//     of << M;
    Test_rcm();
}
catch (DROPSErrCL d) {
    d.handle();
}
    return 0;
}

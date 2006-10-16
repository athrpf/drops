#include "num/spmat.h"
#include "misc/problem.h"
#include <iostream>

using namespace DROPS;


void
depth_first (IdxT v, const MatrixCL& M, PermutationT& p, IdxT& idx);

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

    IdxT idx= 0;
    p.assign( M.num_rows(), NoIdx);
    for (IdxT v= 0; v < M.num_rows(); ++v)
        if (p[v] == NoIdx) depth_first( v, M, p, idx);

    if (idx != M.num_rows())
        throw DROPSErrCL( "topo_sort: Sort failed.\n");

}


void
depth_first (IdxT v, const MatrixCL& M, PermutationT& p, IdxT& idx)
{
    IdxT w;
    --p[v];
    for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
        if (M.val( e) == 0.0 || (w= M.col_ind( e)) == v)
            continue;
        if (p[w] == NoIdx)
            depth_first( w, M, p, idx);
        if (p[w] == NoIdx - 1)
            std::cerr << "Cycle detected at vertex " << w
                << ". Edge (" << v << ", " << w << ") ignored.\n";
    }
    p[v]= idx++;
}


int main ()
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
    std::cerr << "M (the digraph):\n" << M << '\n';

    PermutationT p;
    topo_sort ( M, p);
    seq_out( p.begin(), p.end(), std::cerr);

    VectorCL v( 3);
    v[0]= 0.; v[1]= 1.; v[2]= 2.;
    permute_Vector( v, p);
    std::cerr << v << '\n';

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
    
    std::cerr << "Matrix:\n" << m << '\n' << "Permutation:\n";
    seq_out( p.begin(), p.end(), std::cerr);

    m.permute_columns( p);
    std::cerr << "Columns permuted:\n" << m << '\n';

    m2.permute_rows( p);
    std::cerr << "Rows permuted:\n" << m2 << '\n';

    m3.permute_columns( p);
    m3.permute_rows( p);
    std::cerr << "Columns and rows permuted:\n" << m3 << '\n';

    
    M.permute_rows( p);
    std::cerr << "M, rows permuted:\n" << M << '\n';
    M.permute_columns( p);
    std::cerr << "M permuted:\n" << M << '\n';
    return 0;
}

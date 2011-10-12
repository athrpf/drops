/// \file renumber.tpp
/// \brief implementation of bandwith reduction and downwind numbering (templates)
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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
 * Copyright 2011 LNM/SC RWTH Aachen, Germany
*/

namespace DROPS
{

//=============================================================================
//  Reverse Cuthill-McKee ordering
//=============================================================================

template <typename T>
  void
  breadth_first (size_t v, const SparseMatBaseCL<T>& M, PermutationT& p, size_t& idx,
    std::vector<size_t>& degree)
{
    typedef std::deque<size_t> QueueT;
    QueueT Q;

    // Insert v into queue and number v.
    Q.push_back( v);
    p[v]= idx++;

    std::vector<DegVertT> succ;
    size_t w;
    while (!Q.empty()) {
        v= Q.front();
        Q.pop_front();
        succ.resize( 0);
        succ.reserve( M.row_beg( v + 1) - M.row_beg( v));
        // Find all unnumbered successors of v.
        for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
            if (M.val( e) == 0.0 || (w= M.col_ind( e)) == v)
                continue;
            if (p[w] == NoVert)
                succ.push_back( std::make_pair( degree[w], w));
        }
        // Number successors by increasing degree and enter into queue.
        std::sort( succ.begin(), succ.end(), less1st<DegVertT>());
        for (size_t i= 0; i < succ.size(); ++i) {
            Q.push_back( succ[i].second);
            p[succ[i].second]= idx++;
        }
    }
}

template <typename T>
  size_t
  pseudo_peripheral_node (std::vector<DegVertT>& V, const SparseMatBaseCL<T>& M,
    size_t& max_iter)
{
    PermutationT p_inv( M.num_rows(), NoVert); // Maps vertex v to its index in V.
    for (size_t i= 0; i < V.size(); ++i)
        p_inv[V[i].second]= i;

    // Select a minimal degree vertex.
    size_t v_old= std::min_element( V.begin(), V.end(), less1st<DegVertT>())->second;

    int l_old= 0;
    size_t v, w;

    for (size_t iter= 0; iter < max_iter; ++iter) {
        // Breadth first traversal
        std::vector<int> level( V.size(), -1); // level[p_inv[v]]==-1 means: v is undiscovered.
        typedef std::deque<size_t> QueueT;
        QueueT Q, Qnext;

        // Insert v_old into queue; v_old is in level 0.
        Qnext.push_back( v_old);
        level[p_inv[v_old]]= 0;

        do {
            Q.swap( Qnext);
            Qnext.clear();
            for (QueueT::iterator it= Q.begin(); it != Q.end(); ++it) {
                v= *it;
                // Find all unnumbered successors of v.
                for (size_t e= M.row_beg( v); e < M.row_beg( v + 1); ++e) {
                    if (M.val( e) == 0.0 || (w= M.col_ind( e)) == v)
                        continue;
                    if (level[p_inv[w]] == -1) {
                        Qnext.push_back( w);
                        level[p_inv[w]]= level[p_inv[v]] + 1;
                    }
                }
            }
        }
        while (!Qnext.empty());
        // Now Q contains the last level. It is not empty, if V is not empty.

        // Select a minimal degree vertex from the last level.
        std::vector<DegVertT> lastlevel;
        for (size_t i= 0; i < Q.size(); ++i) {
            lastlevel.push_back( V[p_inv[Q[i]]]);
        }
        v= std::min_element( lastlevel.begin(), lastlevel.end(), less1st<DegVertT>())->second;

        if (level[p_inv[v]] <= l_old) {
            max_iter= iter;
            return v_old;
        }
        l_old= level[p_inv[v]];
        v_old= v;
    }
    return v_old;
}

template <typename T>
void
reverse_cuthill_mckee (const SparseMatBaseCL<T>& M_in, PermutationT& p,
    bool use_indegree)
{
    if (M_in.num_rows() == NoVert)
        throw DROPSErrCL( "reverse_cuthill_mckee: Graph is too large.\n");
    if (M_in.num_rows() != M_in.num_cols())
        throw DROPSErrCL( "reverse_cuthill_mckee: Matrix is not square.\n");

    size_t N= M_in.num_rows();
    SparseMatBaseCL<T> const* M;
    if (use_indegree == true)  M= &M_in;
    else {
        SparseMatBaseCL<T>* M2= new SparseMatBaseCL<T>();
        transpose( M_in, *M2);
        M= M2;
    }

    size_t idx= 0;
    p.assign( N, NoVert);

    std::vector<size_t> degree( N);
    for (size_t r= 0; r < N; ++r)
        degree[r]= M->row_beg( r + 1) - M->row_beg( r);
    std::vector<DegVertT> V;
    V.reserve( N);
    for (size_t i= 0; i < N; ++i)
        V.push_back( std::make_pair( degree[i], i));

    do {
        // Find a start vertex v.
        size_t max_iter= 5;
        size_t v= pseudo_peripheral_node( V, *M, max_iter);
        std::cout << "reverse_cuthill_mckee: p_p_n iterations: " << max_iter << '\n';

        // Traverse the graph starting with v and number the vertices.
        breadth_first( v, *M, p, idx, degree);

        V.clear();
        for (size_t i= 0; i < N; ++i) // Number all components of the graph.
            if (p[i] == NoVert)
                V.push_back( std::make_pair( degree[i], i));
    }
    while (!V.empty());

    if (use_indegree == false) delete M;

    // Revert the order of the vertices.
    for (size_t i= 0; i < N; ++i)
        p[i]= N - 1 - p[i];

    if (idx != N)
        throw DROPSErrCL( "reverse_cuthill_mckee: Could not number all unkowns.\n");
}

//=============================================================================
//  Downwind numbering
//=============================================================================

template <typename T>
  const PermutationT&
  TarjanDownwindCL::number_connected_components (const SparseMatBaseCL<T>& M)
{
    const size_t num_verts= M.num_cols();
    new_index.clear();
    new_index.resize( num_verts, NoVert);
    t_index.clear();
    t_index.resize( num_verts, NoVert);
    low_link.clear();
    low_link.resize( num_verts, NoVert);
    component_.clear();
    component_.resize( num_verts, NoVert);
    component_size_.clear();
    idx= 0;
    tidx= 0;
    stack.clear();
    is_on_stack.clear();
    is_on_stack.resize( num_verts, false);

    for (size_t v= 0; v < num_verts; ++v)
        if (new_index[v] == NoVert)
            number_connected_component_of( v, M);

    return new_index;
}

template <typename T>
  void
  TarjanDownwindCL::number_connected_component_of (size_t v, const SparseMatBaseCL<T>& M)
{
    t_index[v]= low_link[v]= tidx++;
    stack.push_back( v);
    is_on_stack[v]= true;

    const double* succ_val= M.GetFirstVal( v);
    for (const size_t* succ= M.GetFirstCol( v), * rowend= M.GetFirstCol( v + 1); succ != rowend; ++succ, ++succ_val) {
        if (*succ_val <= T()) // Edges in the graph correspond to positive entries. Row i contains the successors of vertex i.
            continue;

        if (t_index[*succ] == NoVert) {
            number_connected_component_of( *succ, M);
            low_link[v]= std::min( low_link[v], low_link[*succ]);
        }
        else if (is_on_stack[*succ])
            low_link[v]= std::min( low_link[v], t_index[*succ]);
    }

    if (low_link[v] == t_index[v]) {
        size_t w;
        component_size_.push_back( 0);
        do {
            w= stack.back();
            component_[w]= component_size_.size() - 1;
            ++component_size_.back();
            is_on_stack[w]= false;
            stack.pop_back();
            new_index[w]= idx++;

        } while (w != v);
    }
}

inline std::vector<size_t> TarjanDownwindCL::component (size_t c) const
{
    std::vector<size_t> ret;
    ret.reserve( component_size()[c]);

    for (size_t i= 0; i < component_.size(); ++i)
        if (component_[i] == c)
            ret.push_back( i);

    return ret;
}

inline void TarjanDownwindCL::stats (std::ostream& os) const
{
    const size_t num_nontrivial= num_components() - std::count( component_size().begin(), component_size().end(), 1);
    os << "# unknowns: " << permutation().size() << '\n'
       << "# components: " << num_components() << " of which " << num_nontrivial << " have more than one unknown.\n"
       << "The largest component has " << *std::max_element(component_size().begin(), component_size().end()) << " unknowns.\n";
}


template <class T>
void
sort_row_entries (SparseMatBaseCL<T>& M)
{
    typedef std::pair<T, size_t> PT;
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

/// \brief Counts the edges within the component 'component' (= positive matrix entries between vertices in 'component') and removes the weakest weak_edge_ratio*100 percent of them.
template <class T>
void remove_weak_edges (SparseMatBaseCL<T>& M, const std::vector<size_t>& component, double weak_edge_ratio)
{
    // Collect all edges within the component
    std::vector<T*> edge;
    for (size_t i= 0; i < component.size(); ++i) {
        T* val= M.GetFirstVal( component[i]);
        for (const size_t* col= M.GetFirstCol( component[i]); col != M.GetFirstCol( component[i] + 1); ++col, ++val)
            if (*val > T() && std::binary_search( component.begin(), component.end(), *col))
                edge.push_back( val);
    }
    // std::cout << edge.size() << " edges within the component.\n";

    // Remove the weakest edges
    const size_t weak_edges= weak_edge_ratio*edge.size();
    typename std::vector<T*>::iterator weak_end= edge.begin() + weak_edges;
    std::nth_element(edge.begin(), weak_end, edge.end(), less_by_ptr<T*>());
    for (typename std::vector<T*>::iterator it= edge.begin(); it != weak_end; ++it)
        **it= T();
}

/// \brief Counts the edges within the component 'component' (= positive matrix entries between vertices in 'component') and removes the weakest weak_edge_ratio*100 percent of them.
/// The routine assumes, that the entries of each row of M are sorted by increasing weight. This is contrary to the linear algebra parts of Drops, which assume an ordering according to increasing column index.
template <class T>
void remove_crosswind_edges (SparseMatBaseCL<T>& M, double crosswind_ratio)
{
    for (size_t i= 0; i < M.num_rows(); ++i) {
        T* val= M.GetFirstVal( i);
        T* valend= M.GetFirstVal( i + 1);
        if (val == valend) // empty row
            continue;
        const double low_limit= *(valend - 1) * crosswind_ratio; // largest value in row times crosswind_ratio
        while (val != valend && *val <= low_limit)
            *val++= T();
    }
}

template <class T>
PermutationT IteratedDownwindCL::downwind_numbering (SparseMatBaseCL<T>& M)
{
    const size_t dim= M.num_rows();

    sort_row_entries( M);
    remove_crosswind_edges( M, crosswind_limit_);

    TarjanDownwindCL re_num;
    bool have_large_components;
    int counter= 0;
    do {
        have_large_components= false;
        re_num.number_connected_components( M);
        // re_num.stats( std::cout);
        for (size_t c= 0; c < re_num.num_components(); ++c)
            if (re_num.component_size()[c] > max_rel_component_size_*dim) {
                // std::cout << "...component " << c << " has " << re_num.component_size()[c] << " vertices.\n";
                // std::cout << "...removing the " << weak_edge_ratio_*100. << " percent weakest edges.\n";
                remove_weak_edges( M, re_num.component( c),  weak_edge_ratio_);
                have_large_components= true;
            }
        ++counter;
    } while (have_large_components);
    re_num.stats( std::cout);
    std::cout << "IteratedDownwindCL::downwind_numbering: " << counter << " iterations.\n";

    return re_num.permutation();
}


} // end of namspace DROPS

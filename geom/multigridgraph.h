/// \file multigridgraph.h
/// \brief classes that constitute the multigridgraph
/// \author LNM RWTH Aachen: Patrick Esser, Nils Gerhard, Joerg Grande; SC RWTH Aachen:

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
#ifndef DROPS_MULTIGRIDGRAPH_H
#define DROPS_MULTIGRIDGRAPH_H 

#include "geom/multigrid.h"
#include <tr1/unordered_map>
#include <set>

namespace DROPS
{
/// \brief class for storage of the colored graph for the OpenMP parallelization 
///
/// class stores a reference to the multigrid an Vectors for parallel iteration over the multigrid
/// it stores the colored graph for different levels if needed
class MultiGridGraphCL{
  private:
    typedef std::vector< std::vector< const TetraCL* > > ColoredGraphT;

    std::map<Uint, std::pair<ColoredGraphT, Uint> > graph_; // level, graph, mgversion
    std::map<Uint, std::pair<ColoredGraphT, Uint> >::iterator graph_iterator; 
    const MultiGridCL& MG_; 
    
  public:

    MultiGridGraphCL( const MultiGridCL & MG);
    ~MultiGridGraphCL(){};
    
    /// \brief Creates colored graph for a level
    void create_graph( Uint level);
    
    /// \brief Checks if there exists a graph for the given level and if its up-to-date, if not it creates a new one
    void update_graph( Uint level);
    
    /// \brief Selects the level an defines an iterator with it 
    void select_graph( Uint level);
    
    /// \brief Returns the number of colors   
    Uint get_num_colors();

    /// \brief Returns the number of tetras in a given color   
    Uint get_num_tetras( Uint i);

    /// \brief Returns a pointer to a tetra for color i and tetra j, it uses the in select_graph defined iterator
    const TetraCL* get_tetra_pointer( Uint i, Uint j);
};

} // end of namespace DROPS

#endif

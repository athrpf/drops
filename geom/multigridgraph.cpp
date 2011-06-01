/// \file multigridgraph.cpp
/// \brief classes that constitute the multigridgraph
/// \author LNM RWTH Aachen: Patrick Esser, Nils Gerhard; SC RWTH Aachen:

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
#include "geom/multigridgraph.h"

namespace DROPS
{

MultiGridGraphCL::MultiGridGraphCL( const MultiGridCL& MG) : MG_( MG) {}
    
void MultiGridGraphCL::create_graph( Uint level){
#ifdef _PAR
    ParTimerCL timer;
#else
    TimerCL timer;
#endif
    timer.Start();
    
    const size_t triangTetraNum = MG_.GetTriangTetra().size( level);
    
    // compute and store all corresponding vertices for every tetra
    std::vector< const TetraCL* > tetraVec( triangTetraNum);
    
    typedef std::set<size_t> tetraNumbT;
    typedef std::tr1::unordered_map< const VertexCL*, tetraNumbT> vertexMapT;
    
    vertexMapT vertexMap( triangTetraNum);
    vertexMapT::iterator JT;
    
    MultiGridCL::const_TriangTetraIteratorCL sit = MG_.GetTriangTetraBegin( level);
    
    // for each tetra the four vertices are stored
    for( size_t k=0; k< triangTetraNum; ++k, ++sit) {
        tetraVec[k] =  &(*sit);
        for( int i=0; i<4; ++i)
            vertexMap[ sit->GetVertex(i) ].insert(k);
    }
    
    // for every tetra compute and store all neighboring tetras in set
    std::vector< tetraNumbT> Tetraliste( tetraVec.size() );
    
    tetraNumbT::iterator tit;
        
    #pragma omp parallel for private( tit)
    for( int j=0; j< (int)tetraVec.size(); j++) //integer loop for openmp
        for( int i=0; i<4; ++i) {
            const tetraNumbT& tmpset = vertexMap.find( tetraVec[j]->GetVertex(i))->second; // list of all tetrahedra which are adjacent to vertex j
            for ( tit = tmpset.begin(); tit != tmpset.end(); ++tit)
                Tetraliste[j].insert( *tit);
        }
        
        vertexMap.clear();
    
    // graph coloring algorithm
    std::vector<int> colorlist( tetraVec.size(), 0); // vector for the color of every tetra
    
    int color_max = 0, new_color = 0;
    
    for( size_t j=0; j< tetraVec.size() ; j++){
        new_color = 0;
        std::set<int> forbidden_colors;
        
        for( tit = Tetraliste[j].begin(); tit != Tetraliste[j].end(); ++tit){
            if( *tit != j)
                forbidden_colors.insert( colorlist[*tit]);
        }
        
        while( forbidden_colors.find( new_color) != forbidden_colors.end())
            new_color++;
        
        if( new_color > color_max)
            color_max =  new_color;
        colorlist[j] = new_color;
    }
    Tetraliste.clear();
    
    ColoredGraphT graph( color_max+1);
    
    for( size_t j=0; j< tetraVec.size(); j++ )
        graph[colorlist[j]].push_back( tetraVec[j]);
    
    colorlist.clear();
    tetraVec.clear();
    
    // detailed list of colored tetras
    //for( size_t h = 0 ; h < Colors.size(); ++h)
    //    std::cout << "Color " << h << " has " << Colors[h].size() << " tetras!" << std::endl;
    //std::cout << std::endl;
    
    
    // tetra sorting for better alignment in cache
    #pragma omp parallel for
    for( int j = 0; j < (int)graph.size(); ++j)
        sort( graph[j].begin(), graph[j].end());
        
    timer.Stop();
    const double duration = timer.GetTime();
        
    std::cout << "creating graph took " << duration << " seconds, use of " << graph.size() << " colors" << std::endl;
    
    graph_[level] = std::make_pair<ColoredGraphT, Uint>( graph, MG_.GetVersion());
}


void MultiGridGraphCL::update_graph( Uint level){
    graph_iterator = graph_.find( level); 
    
    if ( graph_iterator == graph_.end())
        create_graph(level);

    if(( graph_iterator->second).second != MG_.GetVersion()){
        graph_.erase(graph_iterator);
        create_graph(level);
    }           
}

void MultiGridGraphCL::select_graph( Uint level){
    graph_iterator = graph_.find( level);            
}    

Uint MultiGridGraphCL::get_num_colors(){ 
    return ((graph_iterator->second).first).size();
}

Uint MultiGridGraphCL::get_num_tetras( Uint i){ 
    return ((graph_iterator->second).first[i]).size();
}

const TetraCL* MultiGridGraphCL::get_tetra_pointer( Uint i, Uint j){
    return ((graph_iterator->second).first[i])[j];
}


} // end of namespace DROPS
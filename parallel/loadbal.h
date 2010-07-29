/// \file loadbal.h
/// \brief Loadbalancing of tetrahedal multigrids
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier, Timo Henrich

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

#ifndef DROPS_LOADBAL_H
#define DROPS_LOADBAL_H
#include "parallel/parallel.h"
#include "parallel/partime.h"
#include "parallel/parmultigrid.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "parallel/partitioner.h"
#include <map>
#include <set>
#include <iostream>
#include <fstream>

namespace DROPS{

/****************************************************************************
* L O A D  B A L  I T E R A T O R  C L A S S                                *
****************************************************************************/
/// \brief Class for iterating through the multinodes of this proc to set up a
/// reduced, dual graph
///
/// For numbering the tetrahedra for ParMETIS, we have to go through special
/// tetrahedra. This class helps to touch only the right ones
/****************************************************************************
* L O A D  B A L  I T E R A T O R  C L A S S                                *
****************************************************************************/
class LbIteratorCL
{
  public:
    typedef MultiGridCL::TetraIterator IteratorT;     ///< Iteration through the tetraeders of the multigrid

  private:
    MultiGridCL* mg_;                                 // pointer of the multigrid, that is used for iterating
    IteratorT    pos_;                                // The actual position of this iterator
    Uint Level_,                                      // Level actual level
         TriLevel_;                                   // the iterator should only go through this triangulation-level
                                                      // (normaly set to the last level)

  public:
    /// \brief Constructor
    LbIteratorCL (MultiGridCL* mg, IteratorT pos, Uint Level, Uint TriLevel) :
        mg_(mg), pos_(pos), Level_(Level), TriLevel_(TriLevel) {}

    /// \brief Copyconstructor
    LbIteratorCL (const LbIteratorCL &LbI) :
        mg_(LbI.mg_), pos_(LbI.pos_), Level_(LbI.Level_), TriLevel_(LbI.TriLevel_) {}

    // default destructor

    // Test if a tetra is in a loadbalancing set
    inline bool IsInLbSet() const;                    ///< Test, if the tetraeder, on which this iterator points, is in the loadbalancing set
    inline bool IsInLbSet( const TetraCL&) const;     ///< Test if a tetrahedra is in the loadbalancing set see above

    // access of the tetraeder
    TetraCL& operator *  () const { return  *pos_; }  ///< return a reference onto the tetra, this iterator points to
    TetraCL* operator -> () const { return &*pos_; }  ///< return a pointer to a tetra

    // iterate through the elements of the LoadBalance-Set
    inline LbIteratorCL& operator ++ ();              ///< Go to the next multinode tetra (prefix)
    inline LbIteratorCL  operator ++ (int);           ///< Go to the next multinode tetra (suffix)

    // assignement
    LbIteratorCL& operator = (const LbIteratorCL& it) ///< assignment operator
    {
        mg_ = it.mg_; pos_=it.pos_;
        Level_=it.Level_; TriLevel_ = it.TriLevel_;
        return *this;
    }

    // comparing two iterators
    /// \brief Test, if two iterators are the same by comparing the iterators onto the tetras
    friend bool operator == (const LbIteratorCL& It0, const LbIteratorCL& It1)
        { return It0.pos_ == It1.pos_; }
    /// \brief Test, if two iterators are not the same by comparing the iterators onto the tetras
    friend bool operator != (const LbIteratorCL& It0, const LbIteratorCL& It1)
        { return !(It0 == It1); }
};


/****************************************************************************
* L O A D  B A L  C L A S S                                                 *
****************************************************************************/
/// \brief Load balance of a parallel distributed ParMultiGridCL
///
/// This class uses ParMetis to distribute the mulit-grid on several procs
/// by trying to give all procs the same number of tetras and to reduce the
/// edge-cut.
/****************************************************************************
* L O A D  B A L  C L A S S                                                 *
****************************************************************************/
class LoadBalCL
{
  public:
    typedef void (*weight_function)(const TetraCL&, size_t&, const int);
    struct UnkST                                            ///< Storing information about one dof
    { 
        Uint idx;                                           ///< index type, e.g. velocity
        IdxT dof;                                           ///< system number, position in vector
        size_t numUnk;                                      ///< number of unknowns, e.g. 3 for velocity
        UnkST( Uint idx_, IdxT dof_, size_t numUnk_) 
            : idx(idx_), dof(dof_), numUnk(numUnk_) {}
    };

  private:
    typedef std::map<idxtype, Uint> GraphEdgeCT;            // Type that descriebs an edge of the dual, reduced graph: GraphEdgeCT->first: neibhbor, GraphEdgeCT->second: weight
    typedef std::set<idxtype>       GraphNeighborCT;        // Set of neighbors of a node
    typedef std::set<UnkST> UnkWghtListT;                   ///< List of UnkST

    MultiGridCL                   *mg_;                     // Pointer to the multigrid
    PartitionerCL*                partitioner_;             // Pointer to a graph partitioner
    std::vector<const IdxDescCL*> idx_;                     // information about unknowns
    const VecDescCL*              lset_;                    // Eventually use information about interface for loadbalancing
    const BndDataCL<>*            lsetbnd_;                 // Eventually use information about interface for loadbalancing
    Uint static                   TriangLevel_;             // Triangulation level, on which the LoadBalance should be made, normaly set to LastTriangLevel (static for HandlerGather)
    static idxtype*               myfirstVert_;             // first vertex on this proc (static for HandlerGather!)
    static IFT                    FaceIF_;                  // Interface, for compute the adjacencies over Proc-Boundaries
    int                           weightFct_;               // which weighting function should be used

    void InitIF();                                          // Init the Interfaces
    static void CommunicateAdjacency();                     // Communicate the adjacencies over proc boundaries

    LbIteratorCL GetLbTetraBegin( int TriLevel= -1) const;  // iterator to the first tetra in LoadBalSet
    LbIteratorCL GetLbTetraEnd( int TriLevel= -1) const;    // iterator to the last tetra in LoadBalSet
    void CreateNumbering(int TriLevel=-1);                  // Create a global-numbering of the nodes and store the number of numbered tetras in myVerts_

    Uint EstimateAdj();                                     // compute the number of adjacencies on this proc
    std::set<IdxT> UnkOnTetra(const TetraCL&) const;        // get set of unknowns on a tetra
    void AdjUnrefined( TetraCL&, int&, size_t&);            // Set up adjacenzies and weights for a unrefined tetra of LoadBalSet
    void AdjRefined  ( TetraCL&, int&, size_t&);            // Set up adjacenzies and weights for a refined tetra of LoadBalSet

    inline bool CheckForLsetUnk( const TetraCL& t) const;       // Check if level set unknowns are available on a tetrahedron (and its subs)
    void GetWeightRef( const TetraCL&, size_t& wgtpos);         // Weighting function by determining number of children
    void GetWeightLset( const TetraCL&, size_t& wgtpos);        // Weighting function by determining number of intersected subs
    template <typename SimplexT> 
    void UnkOnSimplex(const SimplexT& s, UnkWghtListT& list, const IdxDescCL* idxDesc) const;   // number of dof on a single vertex/edge
    void UnkOnSingleTetra( const TetraCL&, UnkWghtListT& ) const;                               // number of dof on a single tetrahedron
    void GetWeightUnk( const TetraCL&, size_t& wgtpos);         // Weighting function by determining number of degrees of freedom

  public:
    LoadBalCL(MultiGridCL&, int partitioner, int TriLevel=-1, PartMethod meth = KWay, int weightFct=1);// Constructor
    ~LoadBalCL();                                           // Destructor

    void CreateDualRedGraph(bool geom=false);               // Create the dual reduced Graph for ParMETIS on the last triangulation level
    void DeleteGraph();                                     // free all arrays
    PartitionerCL* GetPartitioner() { return partitioner_; }//getter for the partitioner_ object
    void PartitionPar() { partitioner_->PartGraphPar();  }
    void PartitionSer( int master) { partitioner_->PartGraphSer(master); }
    void Migrate();                                         // do the transfers. Call ParMultiGridCL::XferStart(), XferEnd() befor and after calling this procedure
    void RemoveLbNr();                                      // removes all Lb-numbers

    float GetQuality()         const { return partitioner_->GetGraph().ubvec;}  ///< Get quality-parameter
    void  SetQuality(float ub)       { partitioner_->GetGraph().ubvec = ub;}    ///< Set quality-parameter

    void  Append(const IdxDescCL* id) { if (id!=0) idx_.push_back(id); }        ///< Set information about unknown type
    void  RemoveIdx() { idx_.clear(); }                                         ///< Remove information about unknown type
    void  SetWeightFnct( int wfc) { weightFct_=wfc; }                           ///< Set weighting function
    int   GetWeightFnct() const { return weightFct_; }                          ///< Get weighting function

    void SetLset( const VecDescCL& lset, const BndDataCL<>& lsetbnd) { lset_=&lset; lsetbnd_=&lsetbnd;}
    void RemoveLset() { lset_=0; } 

    ///\name for debug purpose
    //{@
    inline int   GetNumAllVerts() const;                    // number of vertices on all procs
    inline int   GetNumLocalAdjacencies() const;            // number of all local adjacencies
    inline Uint  GetEdgeCut() const;                        // Get the edgecut of the last graph partitioning
    inline Uint  GetMovedMultiNodes() const;                // Get the number of all moved multinodes within last migration
    inline float GetImbalanceTolerance() const;             // Get imbalance tolerance for a vertex
    void   ShowGraph(std::ostream&);                        // prints the graph human readable
    friend std::ostream& operator << (std::ostream&,const LoadBalCL&);  // writes onto the ostream in METIS format!
    void   PrintGraphInfo(std::ostream&) const;             // write graph, so that it can be read easily
    //@}

    /// \name handler for DDD as public, so a wrapper can be used
    //{@
    static int HandlerGather( OBJT, void*);              // gather the LoadBalNr of sending proc
    static int HandlerScatter( OBJT, void*);             // scatter the LoadBalNr of revieving proc
    //@}

};

/****************************************************************************
* L O A D  B A L  H A N D L E R  C L A S S                                  *
****************************************************************************/
/// \brief Handler of loadbalancing of a parallel distributed ParMultiGridCL
///
/// Combines different steps to do a load-balancing of a parallel stored
/// MultiGridCL.
/****************************************************************************
* L O A D  B A L  H A N D L E R  C L A S S                                  *
****************************************************************************/
class LoadBalHandlerCL
{
  private:
    LoadBalCL    * lb_;                                 // Load-Balancing-CL
    MultiGridCL  * mg_;                                 // MultiGridCL that should be load balanced
    PartMethod     strategy_;                           // Strategy for computing the graph-partitioning
    bool           xferUnknowns_;                       // flag if unknwons on simplices should be transfered too
    bool           debugMode_;                          // flag if information about each step should be given
    Uint           movedNodes_;                         // number of moved multinodes
    Uint           edgeCut_;                            // number of cutted edges


  public:
    /// \brief Constructor
    LoadBalHandlerCL(MultiGridCL& mg, Partitioner part, float ub=1.00);
    /// \brief Constructor that creates a distributed multigrid
    LoadBalHandlerCL(const MGBuilderCL&, int master, Partitioner part, PartMethod meth=KWay, bool geom=true, bool debug=false);
    /// \brief Destructor, that frees the memory of the LoadBalCL
    ~LoadBalHandlerCL();

    /// \brief Do a complete migration
    void DoMigration();
    /// \brief Distribute a multigrid, that is stored on a single processor
    void DoInitDistribution(int master=Drops_MasterC);

    /// \name Set and get routines
    //@{
    inline PartMethod   GetStrategy()           const;    ///< Returns the strategy that is used for the graph-partitioning
    inline void         SetStrategy(PartMethod);          ///< Set the strategy that is used for the graph-partitioning
    inline bool         GetXferUnknowns()       const;    ///< Returns if unknowns are transfered too within the migration
    inline void         SetXferUnknowns(bool);            ///< Set flag for transfering unknwons (Recursive or Adaptive)
    void SetLset( const VecDescCL& lset, const BndDataCL<>& lsetbnd) { lb_->SetLset(lset, lsetbnd);}    ///< Set a pointer to the level set function, so the number of intersected subs is used for load balancing
    inline bool         GetDebugMode()          const;    ///< Returns if information about each step should be displayed
    inline void         SetDebugMode(bool);               ///< turn off or on the DebugMode
    inline MultiGridCL& GetMG();                          ///< Get a reference on the MultiGridCL
    inline LoadBalCL&   GetLB();                          ///< Get a reference on the LoadBalCL
    inline Uint         GetEdgeCut()            const;    ///< Get number of cutted edges
    inline Uint         GetMovedMultiNodes()    const;    ///< Get number of moved multi nodes
    inline double       GetTetraBalance()       const;    ///< Get ratio between maximal tetra number and minimal tetra number
    //@}
};



// ------------------------------
// I N L I N E  F U N C T I O N S
// ------------------------------


// LbIteratorCL
//=============

bool LbIteratorCL::IsInLbSet() const
{
    return IsInLbSet(*pos_);
}


bool LbIteratorCL::IsInLbSet(const TetraCL& t) const
/// Test, if a tetraeder represents a mutinode of the LoadBalanceSet
{
    if (t.HasGhost() )
        return false;
    if (t.IsUnrefined() )
    {
        if (t.GetLevel()==0)
            return true;
        else
            return false;   // there is a parent, that represents this tetra
    }

    // check for children in this triangulation level, if there is a child, than this tetra is in the LoadBalanceSet
    for (TetraCL::const_ChildPIterator ch(t.GetChildBegin()), chend(t.GetChildEnd()); ch!=chend; ++ch)
        if ((*ch)->IsInTriang(TriLevel_) )
            return true;

    // If non of the above cases happens, this is not in the LoadBalanceSet
    return false;
}


LbIteratorCL& LbIteratorCL::operator ++ ()
/** Increase the position, until we reach an element of the LoadBalanceSet or the
    end of all tetraeders in the triangulation level _TriLevel
*/
{
    do
    {
        ++pos_;                                         // increase pos

        while (pos_ == mg_->GetTetrasEnd( Level_) )     // if we reaches the end of a level
        {
            if (Level_ < TriLevel_)
                pos_= mg_->GetTetrasBegin( ++Level_);   // goto the next level, if in the right triangulation
            else
                return *this;                           // return end
        }
    } while (!IsInLbSet() );        // until reached an element in the LoadBalanceSet
    return *this;                   // return the found element
}


LbIteratorCL LbIteratorCL::operator ++ (int)
{
    LbIteratorCL tmp(*this);
    ++(*this);
    return tmp;
}


// LoadBalCL
//==========

/// \brief Get the number of all nodes, that this class has numbered
int LoadBalCL::GetNumAllVerts() const
{
    Assert( partitioner_->GetGraph().vtxdist,
           DROPSErrCL("LoadBalCL::GetNumAllVerts: Graph has not been set up"),
           DebugLoadBalC);

    return partitioner_->GetGraph().vtxdist[DynamicDataInterfaceCL::InfoProcs()];
}

/// \brief Get the number of local stored adjacencies
int LoadBalCL::GetNumLocalAdjacencies() const{
    return partitioner_->GetGraph().myAdjs;
}

/// \brief Get imbalance tolerance for a vertex
float LoadBalCL::GetImbalanceTolerance() const{
    return partitioner_->GetGraph().ubvec;
}

/// \brief Number of edges that are cut by the last graph partitioning
Uint LoadBalCL::GetEdgeCut() const
{
    return partitioner_->GetGraph().edgecut;
}

/// \brief Number of multinodes that are moved by the last migration
Uint LoadBalCL::GetMovedMultiNodes() const{
    return partitioner_->GetGraph().movedMultiNodes;
}

/// \brief Check if level set unknowns exists on all vertices and edges of a tetrahedron
bool LoadBalCL::CheckForLsetUnk( const TetraCL& t) const
{
    const Uint idx= lset_->RowIdx->GetIdx();
    for ( TetraCL::const_VertexPIterator sit=t.GetVertBegin(); sit!=t.GetVertEnd(); ++sit){
        if ( !(*sit)->Unknowns.Exist(idx))
            return false;
    }
    for ( TetraCL::const_EdgePIterator sit= t.GetEdgesBegin(); sit!=t.GetEdgesEnd(); ++sit){
        if ( !(*sit)->Unknowns.Exist(idx))
            return false;
    }
    return true;
}

/// \brief Put all unknowns on a vertex or edge into a list of dofs
template <typename SimplexT>
void LoadBalCL::UnkOnSimplex( const SimplexT& s, UnkWghtListT& list, const IdxDescCL* idxDesc) const
{
    const IdxT idx   = idxDesc->GetIdx();
    const Uint numUnk= idxDesc->GetNumUnknownsOnSimplex<SimplexT>();
    if ( s.Unknowns.Exist(idx)){
        const IdxT dof= s.Unknowns(idx);
        list.insert( UnkST(idx, dof, numUnk));
        if ( idxDesc->IsExtended() && idxDesc->GetXidx()[dof]!=NoIdx){
//           std::cout << "Using an extended index ("<<idxDesc->GetXidx()[dof]<<") as weight" << std::endl;
           list.insert( UnkST( idx, idxDesc->GetXidx()[dof], numUnk));
        }
    }
}

/// \brief Compare two UnkST
inline bool operator< (const LoadBalCL::UnkST& a, const LoadBalCL::UnkST& b)
{
    if ( a.idx<b.idx)       return true;
    if ( a.dof<b.idx)       return true;
    if ( a.numUnk<b.numUnk) return true;
    return false;
}


// LoadBalHandlerCL
//=================

PartMethod LoadBalHandlerCL::GetStrategy() const{
    return strategy_;
}

void LoadBalHandlerCL::SetStrategy(PartMethod meth)
/// \param[in] meth Methode that is used for computing the graph-partitioning with ParMetis. If meth==Adaptive, then
///            an adaptive recomputation is made, else if meth==Recursive a complete new partitioning is computed
{
    Assert( (meth==Adaptive || meth==Recursive || meth==Identity),
           DROPSErrCL("LoadBalHandlerCL::SetStrategy: No such method known!"),
           DebugLoadBalC);
    strategy_ = meth;
}

bool LoadBalHandlerCL::GetXferUnknowns() const{
    return xferUnknowns_;
}

void LoadBalHandlerCL::SetXferUnknowns(bool xfer){
    xferUnknowns_ = xfer;
}

bool LoadBalHandlerCL::GetDebugMode()const{
    return debugMode_;
}
void LoadBalHandlerCL::SetDebugMode(bool mode){
    debugMode_ = mode;
}

MultiGridCL& LoadBalHandlerCL::GetMG(){
    return *mg_;
}

LoadBalCL& LoadBalHandlerCL::GetLB(){
    return *lb_;
}

Uint LoadBalHandlerCL::GetEdgeCut() const{
    return edgeCut_;
}

Uint LoadBalHandlerCL::GetMovedMultiNodes() const{
    return movedNodes_;
}

double LoadBalHandlerCL::GetTetraBalance() const
{
    size_t mytetras = mg_->GetTetras().size();
    size_t mintetras= ProcCL::GlobalMin(mytetras),
           maxtetras= ProcCL::GlobalMax(mytetras);
    return static_cast<double>(maxtetras)/(static_cast<double>(mintetras));
}

} // end of namespace DROPS

#endif // _LOADBAL_H_

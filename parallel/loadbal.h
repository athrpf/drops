// **************************************************************************
// File:    loadbal.h                                                       *
// Content: Classes to do load balance of a multi-grid                      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen     *
//          Oliver Fortmeier, RZ RWTH Aachen                                *
// Version: 0.1                                                             *
// Date:                                                                    *
// Begin:   Januar, 01st, 2006                                              *
// **************************************************************************
/// \author Oliver Fortmeier
/// \file loadbal.h
/// \brief Loadbalancing of tetrahedal multigrids

#ifndef _LOADBAL_H_
#define _LOADBAL_H_
#include "parallel/parallel.h"
#include <parmetis.h>
#include <metis.h>
// #include "parallel/metispartioner.h"
#include <ddd.h>

#include "parallel/partime.h"
#include "parallel/parmultigrid.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include <map>
#include <set>
#include <iostream>
#include <fstream>

namespace DROPS{

/// \enum PartMethod Tell, which method should be used to compute graph partition problem
enum PartMethod{
    KWay,       ///< Multilevel method
    Recursive,  ///< bisection
    Adaptive,   ///< adaptive recomputation of graph partitioning
    Identity,   ///< Take partition as it is
    NoMig       ///< No migration!
};


/****************************************************************************
* L O A D  B A L  I T E R A T O R  C L A S S                                *
****************************************************************************/
/// \brief Class for iterating through the multinodes of this proc to set up a
/// reduces, dual graph
///
/// For numbering the the tetraes for ParMETIS, we have to go through special
/// tetraeders. This class helps to touch only the right ones
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
    TetraCL& operator *  () const { return  *pos_; }  ///< return a referenze onto the tetra, this iterator points to
    TetraCL* operator -> () const { return &*pos_; }  ///< return a pointer to a tetra

    // iterate through the elements of the LoadBalance-Set
    inline LbIteratorCL& operator ++ ();              ///< Go to the next mulitnode tetra (prefix)
    inline LbIteratorCL  operator ++ (int);           ///< Go to the next mulitnode tetra (suffix)

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
    typedef idxtype* IndexArray;                            ///< tpye of an array of indices

  private:
    typedef std::map<idxtype, Uint> GraphEdgeCT;            // Type that descriebs an edge of the dual, reduced graph: GraphEdgeCT->first: neibhbor, GraphEdgeCT->second: weight
    typedef std::set<idxtype>       GraphNeighborCT;        // Set of neighbors of a node

    MultiGridCL    *mg_;                                    // Pointer to the multigrid
    std::vector<IdxDescCL*> idx_;                           // information about unknowns
    Uint static    TriangLevel_;                            // Triangulationlevel, on which the LoadBalance should be made, normaly set to LastTriangLevel (static for HandlerGather)
    IndexArray     xadj_,                                   // Startingindex in the array _adjncy, where node[i] writes its neighbors
                   adjncy_,                                 // Adjacencies of the nodes
                   vtxdist_,                                // numbers of nodes, that is stored by all procs
                   vwgt_,                                   // weight of the Nodes
                   adjwgt_,                                 // weight of the edges
                   part_;                                   // resulting array, where each node should be send to
    float*         xyz_;                                    // geometric information about the nodes (==barymetric Center of the Tetras)
    float          ubvec_;                                  // quality of graph partitioning
    int            myVerts_;                                // number of vertices on this proc
    int            myAdjs_;                                 // number of Adjacencies
    int            edgecut_;                                // number of edges, that are cut by ParMETIS
    int            movedMultiNodes_;                        // number of multinodes, that are moved by last migration
    static idxtype myfirstVert_;                            // first vertex on this proc (static for HandlerGather!)
    bool           geom_;                                   // flag, that indicates, if the geometrical datas should be used
    static DDD_IF  FaceIF_;                                 // Interface, for compute the adjacencies over Proc-Boundaries
    bool           useUnkInfo_;                             // Use information about unknowns for load balancing


    void InitIF();                                          // Init the Interfaces
    static void CommunicateAdjacency();                     // Communicate the adjacencies over proc boundaries

    LbIteratorCL GetLbTetraBegin( int TriLevel= -1) const;  // iterator to the first tetra in LoadBalSet
    LbIteratorCL GetLbTetraEnd( int TriLevel= -1) const;    // iterator to the last tetra in LoadBalSet
    void CreateNumbering(int TriLevel=-1);                  // Create a global-numbering of the nodes and store the number of numbered tetras in myVerts_

    Uint EstimateAdj();                                     // compute the number of adjacencies on this proc
    std::set<IdxT> UnkOnTetra(const TetraCL&) const;        // get set of unknowns on a tetra
    Uint GetWeight(const TetraCL&) const;                   // get wheight on a tetra
    Uint AdjUnrefined( TetraCL&, int&);                     // Set up adjacenzies and weights for a unrefined tetra of LoadBalSet
    Uint AdjRefined  ( TetraCL&, int&);                     // Set up adjacenzies and weights for a refined tetra of LoadBalSet


  public:
    LoadBalCL(MultiGridCL&, float ub=1.05, int TriLevel=-1);// Constructor
    ~LoadBalCL();                                           // Destructor

    void CreateDualRedGraph(bool geom=false);               // Create the dual reduced Graph for ParMETIS on the last triangulation level
    void DeleteGraph();                                     // free all arrays
    void ParPartKWay();                                     // Computes a partitioning of a parallel distributed graph
    void AdaptRepart(float quality=1000.0);                 // Computes a re-partitioning of a parallel distributed graph
    void SerPartKWay(PartMethod meth=KWay);                 // Compute a graph-partitioning with metis wit specified method
    void IdentityPart();                                    // Set distribution to identity => no balancing
    void Migrate();                                         // do the transfers. Call ParMultiGridCL::XferStart(), XferEnd() befor and after calling this procedure
    void RemoveLbNr();                                      // removes all Lb-numbers

    float GetQuality()         const { return ubvec_; }     ///< Get quality-parameter
    void  SetQuality(float ub)       { ubvec_=ub;}          ///< Set quality-parameter

    void  Append(IdxDescCL* id) { idx_.push_back(id); }     ///< Set information about unknown typ

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
    static int HandlerGather( DDD_OBJ, void*);              // gather the LoadBalNr of sending proc
    static int HandlerScatter( DDD_OBJ, void*);             // scatter the LoadBalNr of revieving proc
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
    LoadBalHandlerCL(MultiGridCL& mg, float ub=1.00);
    /// \brief Constructor that creates a distributed multigrid
    LoadBalHandlerCL(const MGBuilderCL&, int master, PartMethod meth=KWay, bool geom=true, bool debug=false);
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
/** Increase the position, until we reach an element of the LoadBalaceSet or the
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
    Assert(vtxdist_,
           DROPSErrCL("LoadBalCL::GetNumAllVerts: Graph has not been set up"),
           DebugLoadBalC);

    return vtxdist_[DDD_InfoProcs()];
}

/// \brief Get the number of local stored adjacencies
int LoadBalCL::GetNumLocalAdjacencies() const{
    return myAdjs_;
}

/// \brief Get imbalance tolerance for a vertex
float LoadBalCL::GetImbalanceTolerance() const{
    return ubvec_;
}

/// \brief Number of edges that are cut by the last graph partitioning
Uint LoadBalCL::GetEdgeCut() const
{
    return edgecut_;
}

/// \brief Number of multinodes that are moved by the last migration
Uint LoadBalCL::GetMovedMultiNodes() const{
    return movedMultiNodes_;
}

PartMethod LoadBalHandlerCL::GetStrategy() const{
    return strategy_;
}


// LoadBalHandlerCL
//=================

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

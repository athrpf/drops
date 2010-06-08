/// \file fastmarch.h
/// \brief fast marching method for reparametrization
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier
/// \todo Check periodic boundary conditions and reparametrization 

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


#ifndef DROPS_FASTMARCH_H
#define DROPS_FASTMARCH_H

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/spmat.h"
#include "num/interfacePatch.h"
#include "misc/KDtree.h"
#ifdef _PAR
#  include "parallel/interface.h"
#  include "parallel/exchange.h"
#  include "parallel/partime.h"
#endif
#include <set>
#include <memory>

namespace DROPS
{

/// \brief Store all data needed by reparametrization classes
class ReparamDataCL
{
  public:
    enum { Finished= 1, Close= 2, Handled=3, Far=0};    ///< types of vertices (Handled is just for parallel)

  private:
    /// \brief Initialize data structures to handle periodic boundaries
    void InitPerMap();
    /// \brief Determine coordinates of all vertices
    void InitCoord();
    bool gatherPerp;                      ///< Check if perpendicular foots are to be gathered

  public:
    MultiGridCL&             mg;          ///< reference to the multigrid
    VecDescCL&               phi;         ///< reference to the level set function
    VectorCL                 old;         ///< old values of the level set function
    VectorBaseCL<Point3DCL>  coord;       ///< coordinates of all vertices
    VectorBaseCL<byte>       typ;         ///< type of each vertex, i.e., Finished, Close, Far
    VectorBaseCL<Point3DCL*> perpFoot;    ///< perpendicular foots of frontier vertices

    // data for periodic boundaries
    bool                     per;         ///< periodic boundaries are used
    IdxDescCL*               augmIdx;     ///< augmented index for periodic boundaries
    const BndDataCL<>*       bnd;         ///< boundary for periodic data
    VectorBaseCL<IdxT>       map;         ///< mapping of periodic boundary conditions

  public:
    // \brief Allocate memory, store references and init coordinates as well as map periodic boundary dofs
    ReparamDataCL( MultiGridCL& MG, VecDescCL& Phi, bool GatherPerp, bool Periodic=false, const BndDataCL<>* Bnd=0)
        : gatherPerp(GatherPerp), mg( MG), phi( Phi), old( phi.Data),
          coord( Phi.Data.size()), typ( Far, Phi.Data.size()), 
          perpFoot( (Point3DCL*)0, GatherPerp ? Phi.Data.size() : 0),
          per( Periodic), augmIdx( 0), bnd( Bnd), map( 0)
    { InitPerMap(); InitCoord(); }
    /// \brief Delete all perpendicular feet
    ~ReparamDataCL();

    /// \brief for periodic boundary conditions, some mapping is necessary
    inline IdxT Map( IdxT i) const { return i<phi.Data.size() ? i : map[ i-phi.Data.size()]; }
    /// \brief Normalize b onto unit interval [0,1]
    inline void Normalize( double& b) const;
    /// \brief Normalize (b1,b2) onto unit triangle
    inline void Normalize( double& b1, double& b2) const;
    /// \brief Check if perpendicular foots are to be gathered
    inline bool UsePerp() const { return gatherPerp; }
    /// \brief Assign perpendicular foot
    inline void UpdatePerp( const IdxT, const double, const Point3DCL&);
};

inline void ReparamDataCL::Normalize( double& b) const
/// Normalize b onto unit interval [0,1]
{
    if (b<0) b=0;
    else if (b>1) b=1;
}

inline void ReparamDataCL::Normalize( double& b1, double& b2) const
/// Normalize (b1,b2) onto unit triangle.
/// If the barycentric coordinates b1, b2, 1-b1-b2 are non-negative, the point lies within the unit triangle.
/// Otherwise, the nearest point on the triangle's boundary is computed.
{
    const double a=(1-b1-b2)/2., c=b1-b2;
    if (b1>=0 && b2>=0 && a>=0) // point lies within triangle
        return;
    // check whether nearest point is one of the triangle's vertices
    if (b1<=0 && b2<=0) {
        b1=b2=0; return;
    }
    if (b1>=1 && c>=1) {
        b1=1; b2=0; return;
    }
    if (b2>=1 && c<=-1) {
        b1=0; b2=1; return;
    }
    // nearest point lies on one of the triangle's edges
    if (b2<0) { // project on edge b2=0
        b2=0; return;
    }
    if (b1<0) { // project on edge b1=0
        b1=0; return;
    }
    if (a<0) {  // project on edge b2 = 1 - b1
        b1+=a; b2+=a; return;
    }
}

/// \brief Base class to determine distance to the discrete interface \f$\Gamma_\phi\f$ to vertices close to the zero level of level set function
/** Derived classes must provide the function Perform()
    Note that this class do not store any information but only references to data structures
*/
class InitZeroCL
{
  protected:
    ReparamDataCL& data_;                ///< class storing all data
    std::string    name_;                ///< name of the method
    bool           gatherPerp_;          ///< gathering perpendicular foots

  public:
    /// \brief Construct an InitZeroCL
      InitZeroCL( ReparamDataCL& data, const std::string& name)
        : data_( data), name_(name) {}

    /// \brief Destructor deletes numbering for augmented index
    virtual ~InitZeroCL();

    /// \brief Assign each vertex close to the interface the distance
    /** This function must be implemented by each derived class. As result, the type
        of each vertex lying at the zero level must be determined as Finished, and on each
        vertex labeled as Finished  the distance must be set.
     */
    virtual void Perform() = 0;

    /// \brief Get name of the implemented method
    const std::string& GetName() const { return name_; }
};

/// \brief Init zero vertices as finished, do not modify values
class InitZeroNoModCL : public InitZeroCL
{
  public:
    typedef InitZeroCL base;

  public:
    InitZeroNoModCL( ReparamDataCL& data) : base( data, "NoMod") {}

    /// \brief Just mark all vertices at interface as finished
    void Perform();
};

/// \brief Init zero vertices by P1 approximations
/** Initialize the values of vertices lying at the interface. This class
    is not an abstract class with derived classes for projection and scaling
    because, otherwise, many calls to virtual functions would occur.
    \param scale choose between methods, scale==0 projection, scale==1 scaling
 */
template <int scale=0>
class InitZeroP1CL : public InitZeroCL
{
  public:
    typedef InitZeroCL base;

  private:
    VectorCL   sumNormGradPhi_;
    VectorCL   sumVol_;

    /// \brief Compute value of level set function on a child
    inline void ComputeOnChild( IdxT* Numb, int* sign, const ChildDataCL& data, LocalP2CL<>& PhiLoc);

  public:
    InitZeroP1CL( ReparamDataCL& data)
        : base( data, ( scale==1 ? "Scaling" : "P1 projection" )),
          sumNormGradPhi_( scale==1 ? data.phi.Data.size() : 0),
          sumVol_( scale==1 ? data.phi.Data.size() : 0) { }

    void Perform();
};

/// \brief ComputeOnChild is defined for scale==0 or scale==1, otherwise, this standard template is called
template <int scale>
void InitZeroP1CL<scale>::ComputeOnChild( IdxT* Numb, int* sign, const ChildDataCL& data, LocalP2CL<>& PhiLoc)
{
    throw DROPSErrCL("InitZeroP1CL<scale>::ComputeOnChild: scale not defined");
}

/// \brief Computes the distance-function of the level set of ls to the neighboring
///        vertices of the triangulation.
class InitZeroExactCL : public InitZeroCL
{
  public:
    /// \brief Class for computing the distance of a point to a triangle representing the interface
    class DistanceTriangCL
    {
      private:
        Point3DCL       tri_[3];    ///< corners of the triangle
        QRDecompCL<3,2> qr_;        ///< QR-decomposition of edges starting in tri_[0]
        /// \brief Compute distance between line segment and a point
        bool dist (const Point3DCL&, const Point3DCL&, const Point3DCL& p, double&, Point2DCL*);
      public:
        DistanceTriangCL(const Point3DCL tri[3]); ///< Constructor computes QR decompostion
        DistanceTriangCL(const double*);          ///< Copy all data from buffer
        /// \brief Compute distance between p and triangle
        double dist(const Point3DCL& p, byte* locPerp=(byte*)0, Point3DCL* l=(Point3DCL*)(0));
        const Point3DCL& GetTri( size_t i) const { return tri_[i]; }
        const QRDecompCL<3,2>& GetQR() const { return qr_; }
    };

    /// \brief Helper class to determine the position of a DistTriangCL in a vector
    struct DistTriangIndexHelperCL{
        const TetraCL* tetra;       ///< Pointer to intersected tetrahedra
        int            childNum;    ///< number of intersected child
        int            triangNum;   ///< in {0,1}, i.e., intersection may be a triangle or quadrilateral
        DistTriangIndexHelperCL( const TetraCL* t, int ch, int tri)
            : tetra(t), childNum(ch), triangNum(tri) {}
    };

  public:
    typedef InitZeroCL base;
    typedef std::vector<DistanceTriangCL>             DistTriangVecT;
    typedef std::vector<std::vector<size_t> >         DofToDistTriangT;
    typedef std::map<DistTriangIndexHelperCL, size_t> StorePositionMapT;
    typedef std::set<const TetraCL*>                  TetraSetT;
    typedef std::vector<TetraSetT>                    DofToTetraMapT;

  protected:
    DistTriangVecT    distTriang_;      ///< all triangles representing the interface
    StorePositionMapT distTriangPos_;   ///< Get position of a DistanceTriangCL in distTriang_
    DofToTetraMapT    dofToTetra_;      ///< Map an index to a set of intersected tetra
    DofToDistTriangT  dofToDistTriang_; ///< Each dof a set of triangles is associated

    /// \brief Collect all neighbor and neighbor-neighbor tetras of a level set dof
    void InitDofToTetra();
    /// \brief Updates the neighborhood information: traverse all outgoing edges (=tetras) of idx and add them to their end-vertices.
    void insert_neighbor_tetras (IdxT dof, const TetraCL&, Uint ls_idx, const DofToTetraMapT& idx_to_tetra);
    /// \brief Build all DistanceTriangCL
    void BuildDistTriang();
    /// \brief Associate each dof close to the interface the triangles
    void AssociateTriangles();
    /// \brief Determine distances of all dof in the vicinity of the interface
    virtual void DetermineDistances();

    /// \brief Constructor with a name
    InitZeroExactCL( ReparamDataCL& data, const std::string& name)
        : base( data, name) {}

  public:
    InitZeroExactCL( ReparamDataCL& data)
        : base( data, "Exact Distance") {}
    virtual void Perform();
    /// \brief Clean up memory
    virtual void Clean();
    /// \brief Show memory
    virtual void DisplayMem() const;
};


#ifdef _PAR
/// \brief Parallel implementation of InitZeroExactCL
class ParInitZeroExactCL : public InitZeroExactCL
{
  public:
    typedef InitZeroExactCL base;
    typedef std::set<IdxT>  ToSendDistTriangMapT;

  private:
    static ReparamDataCL*         actualData_;          ///< Pointer to actual data, accessible for DDD
    static ToSendDistTriangMapT*  toSend_;              ///< Pointer to toSendDistTriang_, accessible for DDD
    static DofToDistTriangT*      distTriangs_;         ///< Pointer to reference of distance triangles
    static std::map<int, size_t>* actualOffset_;        ///< actual offset
    ToSendDistTriangMapT          toSendDistTriang_;    ///< set of DistanceTriangCL associated to vertices on processor boundaries
    static size_t                 maxTriangsPerDOF_;    ///< number of local distance triangles
    std::map<int, size_t>         offset_;              ///< number of send elements from neighbors

    /// \brief Gather all DistTriangCL for other processes
    void GatherDistTriang();
    /// \brief Communicate triangles
    void CommunicateDistTriang();
    /// \brief Associate new triangles
    void AssociateTrianglesOnProcBnd();

  public:
    ParInitZeroExactCL( ReparamDataCL& data)
        : base( data) {}
    /// \brief Perform all steps determine distance by exact distances on frontier vertices on multiple processes
    void Perform();
    /// \brief Clean up memory
    void Clean();
    /// \brief Show memory
    void DisplayMem() const;

    /// \name Handler for DDD
    //@{
    template<typename SimplexT> static int ExecGatherDistTriang(OBJT);
    template<typename SimplexT> static int HandlerDistTriangGatherPos(OBJT, void*);
    template<typename SimplexT> static int HandlerDistTriangScatterPos(OBJT, void*);
    //@}
};
extern "C" int ExecGatherDistTriangVertexC(OBJT objp);
extern "C" int ExecGatherDistTriangEdgeC(OBJT objp);
extern "C" int HandlerDistTriangGatherPosVertexC(OBJT objp, void* buf);
extern "C" int HandlerDistTriangGatherPosEdgeC(OBJT objp, void* buf);
extern "C" int HandlerDistTriangScatterPosVertexC(OBJT objp, void* buf);
extern "C" int HandlerDistTriangScatterPosEdgeC(OBJT objp, void* buf);
#endif


/// \brief Determine distances of frontier vertices by P2 projection
/// \todo This class is not yet ready to use
class InitZeroP2CL : public InitZeroExactCL
{
  public:
    typedef InitZeroExactCL base;
    typedef InitZeroP2CL    self;
    typedef P2EvalCL< double, BndDataCL<>, VecDescCL> P2EvalT;
    typedef std::vector<std::vector<size_t> >         DofToRepTetraT;
    struct RepTetra{
        LocalP1CL<Point3DCL> G[10];         ///< gradients on tetra
        SMatrixCL<3,3>       T;             ///< transformation to reference tetra
        double               det;           ///< determinate of transformation to reference tetra
        World2BaryCoordCL    w2b;           ///< transform global to local coordinates
        double               valPhi[10];    ///< values of phi on 10 vertices of a tetra
        Point3DCL            coord[10];     ///< coordinates of 10 vertices of a tetra
        Point3DCL            baryCenter;    ///< barycenter of tetra
        inline RepTetra( const TetraCL&, LocalP1CL<Point3DCL>* Gref, const ReparamDataCL&);
    };

protected:
    /// \name Data for minimizing the function f_P2
    //@{
    static RepTetra* actualTetra_;  ///< tetra intersected by interface
    static Uint*     actualVert_;   ///< vertex in tetra
    //@}
    std::vector<RepTetra>            tetras_;           ///< storing all intersected tetras
    std::map<const TetraCL*, size_t> storePosOfTetra_;  ///< mapping a tetra to store position in tetras_
    DofToRepTetraT                   dofToRepTetra_;    ///< map all corresponding tetras to a dof

    /// \brief Build all RepTetra
    void BuildRepTetra();
    /// \brief Associate each dof close to the interface the tetras
    void AssociateTetras();
    /// \brief Determine distances of all dof in the vicinity of the interface
    virtual void DetermineDistances();

public:
    InitZeroP2CL ( ReparamDataCL& data)
        : base( data, "P2") { throw DROPSErrCL("Sorry, InitZeroP2CL not implemented yet!\n"); }
    /// \brief Determine distances of frontier vertices by P2 projection
    virtual void Perform();
    /// \brief Clean up memory
    void Clean();
    /// \brief Show memory
    void DisplayMem() const;

    /// \brief Nonlinear function for optimizer
    static void f_P2(const int *, const double *x, double *fvec, int *);
};


/// \brief Class for propagating the distance from the zero vertices
class PropagateCL
{
  protected:
    ReparamDataCL& data_;               ///< class storing all data
    std::string    name_;               ///< name of the implemented method

  public:
    PropagateCL( ReparamDataCL& data, const std::string& name)
        : data_(data), name_(name) {}
    virtual ~PropagateCL() {}
    /// \brief perform the propagation
    virtual void Perform()=0;
    /// \brief Get name of the implemented method
    const std::string& GetName() const { return name_; }
};

/// \brief Propagate the values by using the fast marching method (FMM)
class FastmarchingCL : public PropagateCL
{
  public:
    typedef std::pair<double, IdxT> DistIdxT;       ///< Helper type for storing close vertices

    /// \brief Class for handling vertices marked as close
    class CloseContCL
    {
      public:
        typedef std::set<DistIdxT> ListT;
        typedef ListT::const_iterator const_iterator;
        typedef ListT::iterator iterator;
      private:
        ListT List_;
      public:
        CloseContCL() { }
        /// \brief returns closest point to interface and deletes this point from List_
        DistIdxT GetNearest() {
            DistIdxT ret = *List_.begin();
            List_.erase( List_.begin());
            return ret;
        }
        // \brief find element
        iterator find(const DistIdxT &a) { return List_.find(a); }
        /// \brief insert
        void insert(const DistIdxT &a) { List_.insert(a); }
        /// \brief remove
        void erase(const DistIdxT &a) { List_.erase( List_.find(a)); }
        /// \brief remove
        void erase(iterator &a) { List_.erase(a); }
        /// \brief Check if list is empty
        bool empty() const { return List_.empty(); }
        /// \brief Get elements in the list
        size_t size() const { return List_.size(); }
    };

    typedef PropagateCL base;                       ///< base class
    typedef SArrayCL<IdxT,4>        ReprTetraT;     ///< representing a tetrahedron
    typedef std::vector<ReprTetraT> VertexNeighT;   ///< neighbor tetrahedra of a DoF

  protected:
    VectorBaseCL<VertexNeighT> neigh_;              ///< neighbor tetrahedra of DoF
    CloseContCL                close_;              ///< vertices marked as close vertices
    size_t                     usedMem_;            ///< Count memory

    /// \brief initialize neighbor vertices
    virtual void InitNeigh();
    /// \brief initialize close set
    virtual void InitClose();
    /// \brief Update value on a vertex and put this vertex into close set
    void Update( const IdxT);
    /// \brief Compute projection on linearized level set function on child
    double CompValueProj( IdxT Nr, int num, const IdxT upd[3]) const;
    /// \brief Compute the distances
    void DetermineDistances();

  public:
    FastmarchingCL( ReparamDataCL& data)
        : base( data, "Fast-Marching-Method") {}
    virtual ~FastmarchingCL() {}
    /// \brief Determine unsigned distances by the fast marching method
    virtual void Perform();
};

#ifdef _PAR
/// \brief Performing the FMM on a master process
class FastmarchingOnMasterCL : public FastmarchingCL
{
  public:
    typedef FastmarchingCL base;
    /// \brief structure for sending mark and value
    struct CoupMarkValST{
        byte mark;      ///< mark of dof
        double val;     ///< value of level set function
    };

  private:
    static ReparamDataCL*      actualData_; ///< make data accessible for DDD

    static std::vector<IdxT>   globNumb_;   ///< mapping local DoF on global DoF number
    static std::vector<IdxT>   locNumb_;    ///< mapping global DoF number on local exclusive DoF number
    std::vector<int>           offset_;     ///< offset for each process
    std::vector<int>           exclusive_;  ///< number of exclusive dof per process
    size_t                     size_;       ///< number of level set dof
    int                        master_;     ///< process that performs the FMM

    /// \brief Create global numbering for fast marching algorithm
    void CreateGlobNumb();
    /// \brief Communicate Finished marked DoF
    void DistributeFinished();
    /// \brief Collect all local data
    void CollectLocalData(std::vector<byte>& typ, std::vector<double>& coord, std::vector<double>& values, std::vector<IdxT>& tetraList) const;
    /// \brief Send all information to master-process
    void Collect();
    /// \brief Master send all information to other processes
    void Distribute();
    /// \brief Clean up memory
    void CleanUp();

  public:
    FastmarchingOnMasterCL( ReparamDataCL& data, int master=Drops_MasterC)
        : base( data), master_(master) {}
    ~FastmarchingOnMasterCL() { CleanUp(); }
    /// \brief Send data to master and perform the FMM on this process
    void Perform();

  public:
    /// \brief Gather global number of DoFs
    template<class SimplexT> static int HandlerGlobDOFGather(OBJT, void*);
    /// \brief Scatter global number of DoFs
    template<class SimplexT> static int HandlerGlobDOFScatter(OBJT, void*);
    /// \brief Gather finished DoF
    template<class SimplexT> static int HandlerFinishedGather(OBJT, void*);
    /// \brief Scatter finished DoF
    template<class SimplexT> static int HandlerFinishedScatter(OBJT, void*);
};
/// \name  Wrapper for DDD
//@{
extern "C" int HandlerGlobDOFGatherVertexC(OBJT objp, void* buf);
extern "C" int HandlerGlobDOFGatherEdgeC(OBJT objp, void* buf);
extern "C" int HandlerGlobDOFScatterVertexC(OBJT objp, void* buf);
extern "C" int HandlerGlobDOFScatterEdgeC(OBJT objp, void* buf);
extern "C" int HandlerFinishedGatherVertexC(OBJT objp, void* buf);
extern "C" int HandlerFinishedGatherEdgeC(OBJT objp, void* buf);
extern "C" int HandlerFinishedScatterVertexC(OBJT objp, void* buf);
extern "C" int HandlerFinishedScatterEdgeC(OBJT objp, void* buf);
//@}
#endif

/// \brief Determine distances by direct distance computing to frontier vertices and
///        perpendicular feet
class DirectDistanceCL : public PropagateCL
{
  public:
    typedef PropagateCL base;

  protected:
    KDTreeCL<double>* kdTree_;      ///< k-d tree to search for nearest neighbors
    size_t            numNeigh_;    ///< number of neighbors
    VectorCL          front_;       ///< coordinates of frontier vertices and perpendicular feet
    VectorCL          vals_;        ///< values of phi on frontier vertices and perpendicular feet

    /// \brief Initialize the vectors front_ and vals_
    void InitFrontVector();
    /// \brief Build the KD-Tree of front_
    void BuildKDTree();
    /// \brief Determine distances by using the kd-tree
    void DetermineDistances();

  public:
    DirectDistanceCL( ReparamDataCL& data, size_t numNeigh=100)
        : base( data, "Direct Distance"), kdTree_(0), numNeigh_( numNeigh), front_(0), vals_(0) {}
    /// \brief Determine unsigned distances by the direct computing of distances
    virtual void Perform();
    /// \brief Show memory
    void DisplayMem() const;
};

#ifdef _PAR
/// \brief Determine distances by direct distance computing to frontier vertices and
///        perpendicular feet on multiple processes
class ParDirectDistanceCL : public DirectDistanceCL
{
  public:
    typedef DirectDistanceCL base;
    struct TransferST{ double value; Point3DCL perp; int procID; };
    /// \brief mapping of frontier vertices to a list of values on all
    ///        processes storing the frontier vertex
    typedef std::map<IdxT, std::list<TransferST> > MultiFrontT;

  private:
    static MultiFrontT    onProc_;      ///< static, so DDD may access container
    static ReparamDataCL* actualData_;  ///< static, so DDD may access data

    /// \brief Communicate values and perpendicular feet on frontier vertices located at
    ///    process boundaries
    void CommunicateFrontierSetOnProcBnd();
    /// \brief Communicate all values and perpendicular feet of all processes
    void GatherFrontier();
    /// \brief Clean up memory
    void CleanUp();

  public:
    ParDirectDistanceCL( ReparamDataCL& data, size_t numNeigh=100)
        : base(data, numNeigh) {}
    ~ParDirectDistanceCL() {CleanUp(); }
    /// \brief Determine unsigned distances by the direct computing of distances on multiple processors
    void Perform();

  public:
    /// \name Handlers for DDD
    //@{
    template<class SimplexT>
      static int HandlerFrontierGather(OBJT, void*);         ///< Gather finished DoF
    template<class SimplexT>
      static int HandlerFrontierScatter(OBJT, void*);        ///< Scatter finished DoF
    //@}
};
/// \name  Wrapper for DDD
//@{
extern "C" int HandlerFrontierGatherVertexC(OBJT objp, void* buf);
extern "C" int HandlerFrontierGatherEdgeC(OBJT objp, void* buf);
extern "C" int HandlerFrontierScatterVertexC(OBJT objp, void* buf);
extern "C" int HandlerFrontierScatterEdgeC(OBJT objp, void* buf);
//@}
#endif

// fwd declaration
class ReparamFactoryCL;

/// \brief Perform the reparametrization
class ReparamCL
{
  public:
    friend class ReparamFactoryCL;

  private:
    ReparamDataCL data_;            ///< data for reparametrization
    InitZeroCL*   initZero_;        ///< initialization of frontier vertices
    PropagateCL*  propagate_;       ///< propagate values of the level set function

    /// \brief Restore all signs
    void RestoreSigns();

  public:
    /// \brief Constructor
    ReparamCL( MultiGridCL& mg, VecDescCL& phi, bool gatherPerp, bool periodic=false, const BndDataCL<>* bnd=0);
    ~ReparamCL();
    /// \brief Perform the reparametrization
    void Perform();
};

/// \brief Generating a reparametrization class
/** Chose by parameter method one of the following methods
    <table border="3">
    <tr><th> method </th><th> Init Zero Method  </th><th> Propagation method            </th></tr>
    <tr><td>  00    </td><td> No modification   </td><td> Fast marching method          </td></tr>
    <tr><td>  01    </td><td> P1 Scaling        </td><td> Fast marching method          </td></tr>
    <tr><td>  02    </td><td> P1 projection     </td><td> Fast marching method          </td></tr>
    <tr><td>  03    </td><td> Exact Distance    </td><td> Fast marching method          </td></tr>
    <tr><td>  10    </td><td> No modification   </td><td> Direct distance with KD trees </td></tr>
    <tr><td>  11    </td><td> P1 Scaling        </td><td> Direct distance with KD trees </td></tr>
    <tr><td>  12    </td><td> P1 projection     </td><td> Direct distance with KD trees </td></tr>
    <tr><td>  13    </td><td> Exact Distance    </td><td> Direct distance with KD trees </td></tr>
    </table>
*/
class ReparamFactoryCL
{
  public:
    ReparamFactoryCL() {}
    /// \brief Construct a reparametrization class
    static std::auto_ptr<ReparamCL> GetReparam( MultiGridCL& mg, VecDescCL& phi, int method=03, bool periodic=false, const BndDataCL<>* bnd=0);
};

} // end of namespace DROPS

#include "levelset/fastmarch.tpp"

#endif

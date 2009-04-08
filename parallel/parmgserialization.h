//**************************************************************************
// File:    parmgserializationcl.h                                         *
// Content: Definition of Class ParMGSerializationCL                       *
// Author:  Oliver Fortmeier, SC RWTH Aachen                               *
//          Timo Henrich, SC RWTH Aachen
// Version: 0.1                                                            *
// History: begin - November, 12 2007                                      *
//**************************************************************************
#ifndef DROPS_PARSERIAL_H
#define DROPS_PARSERIAL_H

#include "geom/multigrid.h"
#include "num/bndData.h"
#include "geom/topo.h"
#include "misc/utils.h"
#include "misc/problem.h"
#include <istream>
#include <vector>
#include <map>
#include <cmath>

namespace DROPS
{

/*******************************************************************
*   P A R  S E R I A L I Z E  E R R  C L                           *
*******************************************************************/

/// \brief Error handling of parallel serialization
class ParSerializeErrCL : public DROPSErrCL
{
  public:
    typedef DROPSErrCL base;
    using base::_ErrMesg;

  private:
    int level_;                 // Level of error (0 harmless, 1 serious)

  public:
    ParSerializeErrCL(const std::string& str, const int& lvl)
        : base(str), level_(lvl) {}
    ParSerializeErrCL(const ParSerializeErrCL& err)
        : base(err._ErrMesg), level_(err.level_) {}
    ~ParSerializeErrCL() { }

    /// \brief Handle errors, that arises while serialize or read data
    void handle() const{
        if (level_>0)
            base::handle();
        else{
            base::what(std::cout);
            if (ProcCL::IamMaster())
                std::cout << "===> Skipping serialization and continue program <===" << std::endl;
        }
    }
};


/*******************************************************************
*   P A R M G S E R I A L I Z A T I O N   C L                      *
*******************************************************************/

/// \brief Parallel-version of class MGSerializationCL
class ParMGSerializationCL
/** Every worker collects its edges,vertices,faces,tetras etc. from the given
    multigrid and sends the data back to the master process. The master process
    itself collects the data from all proces and writes it into seperate files.
    <p>
    Therefore four classes for storing information about all types of simplices
    are used, i.e. VertexInfoCL, EdgeInfoCL, FaceInfoCL and TetraInfoCL which
    are derived by the base class BasicGeomCL. To transfer and handle the
    boundary- as well as the children-information two more classes are used,
    i.e. BndVertexInfoCL and TetraChildsInfoCL.
    <p>
    The worker processes collect only information about vertices, edges and
    tetrahedra, if they own an exclusive copy of the simplex. Due to the fact,
    that neighbor information of tetrahedra are stored by the faces, all faces
    are transmitted to the master process, which put neighborhood information
    to the right places within the faces.
    <p>
    This class is also capable to write out degrees of freedom (dof) on vertices
    and edges. This can be done by the function \a WriteDOF. This procedure
    follows the same scheme as the procedures to serialize the multigrid.
*/
{
  public:
    /// \brief Basic storage class for all simplices
    class BasicGeomCL
    /** All derived classes are used to send information from the worker
        processes to the master one*/
    {
      public:
        Uint gid;       ///< global id of simplex
        Uint level;     ///< level of the simplex

      public:
        /// \brief Predicate used for sorting elements in right order.
        static bool compareLevel( const BasicGeomCL &a, const BasicGeomCL &b ) {
            return a.level < b.level;
        }
    };

    /// \brief Storage class for vertices
    class VertexInfoCL : public BasicGeomCL
    {
      public:
        typedef BasicGeomCL base;
        using base::gid;
        using base::level;

      public:
        double choords[3];          ///< Coordinates
        int    remove;              ///< Remove mark
        int    isOnBoundary;        ///< flag if vertex lies on boundary
        Usint  bndIdxStart;         ///< first boundary segment
        Usint  bndIdxEnd;           ///< last boundary segment

      public:
        /// \brief Default Constructor
        VertexInfoCL(){}
        /// \brief Constructor with a vertex
        VertexInfoCL(const VertexCL &src){
            Init(src);
        }

        /// \brief Initialitzes this object with a given vertex.
        void Init(const VertexCL&);
    } ;

    /// \brief Storage class for edges
    class EdgeInfoCL : public BasicGeomCL
    {
      public:
        typedef BasicGeomCL base;
        using base::gid;
        using base::level;

      public:
        Uint vertices[2];           ///< global ids of the both vertices
        Uint midVertex;             ///< global id of the midvertex
        unsigned short int accMFR;  ///< accumulaed "mark for refinement"
        Usint bndIdxStart;          ///< first boundary segment index
        Usint bndIdxEnd;            ///< second boundary segment index
        Uint remove;                ///< Remove mark

      public:
        /// \brief Default-Constructor
        EdgeInfoCL(){}

        /// \brief Constructor with an edge
        EdgeInfoCL(const EdgeCL &src){
            Init(src);
        }

        /// \brief Initialitzes this object with a given edge.
        void Init(const EdgeCL&);
    };

     /// \brief Storage class for faces
    class FaceInfoCL : public BasicGeomCL
    {
      public:
        typedef BasicGeomCL base;
        using base::gid;
        using base::level;

      public:
        Uint neighbor[4];           ///< neighbors (up to four at green ref rules)
        unsigned short int bndIdx;  ///< boundary index, if face lies on a boundary
        Uint remove;                ///< Remove mark
        Uint isExclusive ;          ///< flag if this face has to be written out

      public:
        /// \brief Constructor
        FaceInfoCL(){}

        /// \brief Constructor with a face
        FaceInfoCL(const FaceCL &src){
            Init(src);
        }

        /// \brief Initialitzes this object with a given face
        void Init(const FaceCL&);
    };

    class TetraInfoCL : public BasicGeomCL
    {
      public:
        typedef BasicGeomCL base;
        using base::gid;
        using base::level;

      public:
        Uint vertices[4];   ///< four vertices of a tetrahedron
        Uint edges[6];      ///< six edges of a tetrahedron
        Uint faces[4];      ///< four faces of a tetrahedron
        Uint parentTetra;   ///< parent tetrahedron, if tetrahedron is refined
        Uint refRule;       ///< refinement rule
        Uint refMark;       ///< refinement mark

      public:
        /// \brief Constructor
        TetraInfoCL(){}

        /// \brief Constructor with a tetrahedron
        TetraInfoCL(const TetraCL &src){
            Init(src);
        }

        /// \brief Initialitzes this object with a given tetra
        void Init(const TetraCL&);
    };

    /// \brief Storage class for boundary vertices
    class BndVertexInfoCL
    {
      public:
        Uint gid;                   ///< global id of the vertex
        Uint bndIdx;                ///< index of the boundary segment
        double choords[2];          ///< 2D Coordinates within the boundary segment

      public:
        /// \brief Constructor
        BndVertexInfoCL(){}

        /// \brief Constructor with a vertex and boundary information
        BndVertexInfoCL(const VertexCL &srcVertex, const BndPointCL &srcBndVertex){
            Init(srcVertex,srcBndVertex);
        }

        /// \brief Initialitzes this object with a given vertex and vertex-boundary
        void Init(const VertexCL&, const BndPointCL&);
    };

    /// \brief Storage class for child information
    class TetraChildsInfoCL
    {
      public:
        Uint gid;                   ///< global id of parent tetrahedron
        Uint childs[8];             ///< global ids of all children

      public:
        /// \brief Constructor
        TetraChildsInfoCL() {}

        /// \brief Constructor with a given tetrahedron
        TetraChildsInfoCL(const TetraCL &src){
            Init(src);
        }

        /// \brief initialitzes this object with a given tetrahedron
        void Init(const TetraCL&);
    };

    /// \brief Storage class for a scalar dof value.
    class DofScalarCL
    {
        public:
            Uint gid;
            double data;
    };

    /// \brief Storage class for a vectorial dof value.
    class DofVecCL
    {
        public:
            Uint gid;
            double data[3];
    };

  protected:
    // Path and/or File-Prefix
    std::string  path_;

  private:
    // Reference to the used multigrid
    const MultiGridCL& mg_;

    // All used MPI datatype
    ProcCL::DatatypeT vertexStMPI_;
    ProcCL::DatatypeT bndVertexStMPI_;
    ProcCL::DatatypeT edgeStMPI_;
    ProcCL::DatatypeT faceStMPI_;
    ProcCL::DatatypeT tetraStMPI_;
    ProcCL::DatatypeT tetraChildsStMPI_;
    ProcCL::DatatypeT dofScalarMPI_;
    ProcCL::DatatypeT dofVecMPI_;

    // Buffer for all known elements
    std::vector<VertexInfoCL> vertexBuffer_;
    std::vector<BndVertexInfoCL> bndVertexBuffer_;
    std::map<Uint,std::vector<BndVertexInfoCL> > bndVertexMap_;
    std::vector<EdgeInfoCL>   edgeBuffer_;
    std::vector<FaceInfoCL>   faceBuffer_;
    std::vector<TetraInfoCL>  tetraBuffer_;
    std::vector<TetraChildsInfoCL>  tetraChildsBuffer_;

    // Buffers for local Dofs
    std::vector<DofScalarCL>    dofScalBuffer_;
    std::vector<DofVecCL>       dofVecBuffer_;


    std::map<Uint,double>       dofScalMap_;
    std::map<Uint,double*>      dofVecMap_;

    // The adress map
    // maps an new gid to the next not used integer
    typedef std::map<Uint,Uint> MapIdToGIDT;
    MapIdToGIDT addrMapVertex_;
    MapIdToGIDT addrMapEdge_;
    MapIdToGIDT addrMapFace_;
    MapIdToGIDT addrMapTetra_;

    // Tags to identify messages
    const int vertexTag_    ;
    const int bndVertexTag_ ;
    const int edgeTag_      ;
    const int faceTag_      ;
    const int tetraTag_     ;
    const int tetraChildTag_;
    const int dofScalarTag_;
    const int dofVecTag_;

    // Id of the master proc
    int masterProc_;

    // Check if an edge or vertex is owned by proc with lowest rank
    template<typename SimplexT>
    bool AmIResponsibleProc(const SimplexT& s) { return s.IsExclusive(PrioNeutral);}

    // Initializes the pAddrMap. The function is called in WriteMG by the master-process.
    void PrepareData();
    // helper function for PrepareData in order to register faces
    void RegisterFace(FaceInfoCL&);
    // helper function for PrepareData in order to store parent on child in neighbors on the same side
    void MakeNeighborsConsistent(FaceInfoCL&);

    // Get internal id of a global id
    Uint AddrMap_getLocalId(MapIdToGIDT &tmpMap,Uint searchGID);

    // Map gid to a consecutive internal id
    bool AddrMap_storeGID(MapIdToGIDT &tmpMap,Uint newGID);


    // Routines for sending and moving objects back to the
    // master-proc.
    void SendEdges    ();
    void SendFaces    ();
    void SendVertices ();
    void SendTetras   ();
    void MoveEdges    ();
    void MoveFaces    ();
    void MoveVertices ();
    void MoveTetras   ();

    // Write onjects into files
	void WriteVertices(const std::string &path);
    void WriteEdges   (const std::string &path);
    void WriteFaces   (const std::string &path);
    void WriteTetras  (const std::string &path);

   // Receiving-Routine
   void FetchData();

   // Register-Functions for new MPI-datatypes
   void CreateVertexMPI();
   void CreateBndVertexMPI();
   void CreateEdgeMPI();
   void CreateFaceMPI();
   void CreateTetraMPI();
   void CreateTetraChildsMPI();
   void CreateDofScalarMPI();
   void CreateDofVecMPI();

   // Collect information about dofs on vertices and edges.
   void CollectDOF(const VecDescCL*);
   // Send dofs to master processor
   void SendDOF();
   // Free used local buffers
   void FreeLocalDOFBuffer();

   // Recive dof on master processor
   void RecieveDOF(const VecDescCL* vec, int p);

    /**
    * \brief Expand the current DOF-Map with dof supplied in buffer.
    * \param pBuffer A buffer containing scalar dof values
    * \param bufferSize Number of elements in pBuffer
    **/
    void ExpandDOFMap(const DofScalarCL* pBuffer,Uint bufferSize );

    /**
    * \brief Expand the current DOF-Map with dof supplied in buffer.
    * \param pBuffer A buffer containing vectorial dof values
    * \param bufferSize Number of elements in pBuffer
    **/
	void ExpandDOFMap(const DofVecCL* pBuffer,Uint bufferSize );

  public:
  	/// \brief Constructor
    ParMGSerializationCL(const MultiGridCL& mg, const std::string& path, int master=Drops_MasterC)
        : path_(path), mg_(mg),
          vertexStMPI_(ProcCL::NullDataType), bndVertexStMPI_(ProcCL::NullDataType),
          edgeStMPI_(ProcCL::NullDataType), faceStMPI_(ProcCL::NullDataType),
          tetraStMPI_(ProcCL::NullDataType), tetraChildsStMPI_(ProcCL::NullDataType),
          dofScalarMPI_(ProcCL::NullDataType), dofVecMPI_(ProcCL::NullDataType),
          dofScalBuffer_(), dofVecBuffer_(), dofScalMap_(), dofVecMap_(),
          vertexTag_(100),bndVertexTag_(110), edgeTag_(200),
          faceTag_(300), tetraTag_(400), tetraChildTag_(410),
          dofScalarTag_(420),dofVecTag_(430), masterProc_(master)
    {
        // Register all needed mpi-types
        CreateVertexMPI();
        CreateBndVertexMPI();
        CreateEdgeMPI();
        CreateFaceMPI();
        CreateTetraMPI();
        CreateTetraChildsMPI();
    }
    ~ParMGSerializationCL();

    /// \brief Frees the memory of all used buffers
    void Clear();

    /// \brief Print information about written objects
    void printStat();

    /// \brief Write MultiGridCL into files, so a MGBuilderCL can use them
    /// \param path An alternative output-path
    void WriteMG (const std::string path="");

    /// \brief Write dofs into a file
    void WriteDOF(const VecDescCL*, const std::string&,const std::string path ="");
};

/// \brief Read dof from a file
void ReadDOF(const MultiGridCL&, VecDescCL*, const std::string&);


/// \brief Class for write out the solution of a two-phase flow problem
template <typename StokesT, typename LevelsetT>
class TwoPhaseSerializationCL : public ParMGSerializationCL
{
  public:
    typedef ParMGSerializationCL                        base;
    typedef TwoPhaseSerializationCL<StokesT, LevelsetT> self;

  private:
    bool        overwrite_;
    StokesT&    Stokes_;
    int         localTimeStep_;
    std::string format_;
    int         digits_;
    LevelsetT&  Levelset_;
    using base::path_;

  public:
    TwoPhaseSerializationCL(const MultiGridCL& mg, const std::string& path,
                            StokesT& Stokes, LevelsetT& Levelset,
                            bool overwrite=true, Uint maxtimesteps=100)
      : base(mg, path, ProcCL::Master()), overwrite_(overwrite), Stokes_(Stokes), localTimeStep_(0),
        format_("%0*d"),
        digits_((int)std::ceil(std::log((double)(maxtimesteps!=0 ? maxtimesteps+1 : 2))/std::log((double)10))),
        Levelset_(Levelset) {}

    /// \brief Write information about geometry, pressure, velocity and level-set into files
    void Serialize(int timestep=-1);
};


// *****************************************************************************
// DECLARATION OF TEMPLATE AND INLINE FUNCTIONS
// *****************************************************************************

template <typename StokesT, typename LevelsetT>
  void TwoPhaseSerializationCL<StokesT, LevelsetT>::Serialize(int timestep)
/** This function writes all information about unknowns and geometry into files
    in a single directory. If the flag overwrite is not set, this function
    creates for each serialization a new directory "TimeStep_xxx".
    \param timestep default -1: own counter is used to create unique directories.
    \pre path_ (set in the base class) must be an existing directory
*/
{
    try{
        std::string old_path=path_;
        if (!overwrite_){
            char buffer[72];
            std::sprintf(buffer, format_.c_str(), digits_, (timestep<0 ? localTimeStep_++ : timestep));
            std::string new_path= path_+"/TimeStep_"+buffer+"/";
            if (ProcCL::IamMaster()){
                if (CreateDirectory(new_path)) throw ParSerializeErrCL("TwoPhaseSerializationCL::Serialize: Cannot create directory", 0);
                path_=new_path;
            }
        }
        base::WriteMG();
        base::WriteDOF(&(Stokes_.v),     "velocity");
        base::WriteDOF(&(Stokes_.p),     "pressure");
        base::WriteDOF(&(Levelset_.Phi), "level-set");
        base::Clear();
        path_= old_path;
    }
    catch (ParSerializeErrCL err){
        err.handle();
    }
}

} // End Namespace Drops
#endif

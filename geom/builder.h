/// \file builder.h
/// \brief MGBuilderCL objects for some domains
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Peters, Volker Reichelt; SC RWTH Aachen:

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

/// Remarks: We should use the const-qualifier to make it difficult to
///          accidentally change the multigrid structure from anywhere
///          outside of the multigrid algorithms.
///          Thus the pointer to user data structures should probably be
// /         a pointer to mutable.

#ifndef DROPS_BUILDER_H
#define DROPS_BUILDER_H


#include "geom/multigrid.h"
#include "num/bndData.h"
#include "misc/utils.h"
#include <istream>
#include <map>

#ifdef _PAR
#include "parallel/pardistributeddata.h"
#endif

namespace DROPS
{

/// \brief Class for building a brick
class BrickBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    // for vertices:
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    // for tetras:
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

  protected:
    /// \brief build boundary
    void buildBoundary(MultiGridCL* mgp) const;

  public:
    BrickBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&, Uint, Uint, Uint);

    /// \brief Build a triangulation of a brick
    virtual void
    build(MultiGridCL*) const;
};

/// \brief Class for building a brick with a hole
class CavityBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    const SArrayCL<Uint, 3> cavityorigin_;
    const SArrayCL<Uint, 3> cavity_;
    // for vertices:
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    // for tetras:
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

  protected:
    /// \brief build boundary
    void buildBoundary(MultiGridCL* mgp) const;

  public:
    CavityBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&, Uint, Uint, Uint, SArrayCL<Uint, 3> cavityorigin, SArrayCL<Uint, 3> cavity);

    /// \brief Build a triangulation of a brick
    virtual void
    build(MultiGridCL*) const;
};


/// \brief Class for setting up data structures of a multigrid of a brick shaped
///        domain
/** In order to create a brick on a parallel computer, the following strategy is
    used. Only the master processor creates the multigrid and all other
    processors have to create an "empty" multigrid, i.e. they only create the
    level views and the boundary. Therefore this class is instantiated from the
    non-master processors. The master process instantiates the base class of this
    class.
*/
class EmptyBrickBuilderCL : public BrickBuilderCL
{
  private:
    typedef BrickBuilderCL base_;
    Uint _numLevel;

  public:
    EmptyBrickBuilderCL( const Point3DCL& orig, const Point3DCL& e1,
                         const Point3DCL& e2, const Point3DCL& e3, Uint Level=1)
      : base_( orig, e1, e2, e3, 0, 0, 0), _numLevel( Level) {}

    /// \brief Build level views and boundary information
    virtual void build( MultiGridCL*) const;
};

class EmptyCavityBuilderCL : public CavityBuilderCL
{
  private:
    typedef CavityBuilderCL base_;
    Uint _numLevel;

  public:
    EmptyCavityBuilderCL( const Point3DCL& orig, const Point3DCL& e1,
                         const Point3DCL& e2, const Point3DCL& e3, Uint Level, SArrayCL<Uint, 3>, SArrayCL<Uint, 3>)
      : base_( orig, e1, e2, e3, 0, 0, 0, SArrayCL<Uint, 3>(0u), SArrayCL<Uint, 3>(0u)), _numLevel( Level) {}

    /// \brief Build level views and boundary information
    virtual void build( MultiGridCL*) const;
};

class LBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    const Uint _b1;
    const Uint _b2;
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

 protected:
    void buildBoundary(MultiGridCL* mgp) const;

  public:
    LBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&,
               Uint, Uint, Uint, Uint, Uint);

    virtual void
    build(MultiGridCL*) const;
};


class BBuilderCL : public MGBuilderCL
{
  private:
    const Point3DCL _orig;
    const Point3DCL _e1;
    const Point3DCL _e2;
    const Point3DCL _e3;
    const Uint _n1;
    const Uint _n2;
    const Uint _n3;
    const Uint _b1;
    const Uint _b2;
    const Uint _b3;
    Uint v_idx(Uint i, Uint j, Uint k)         const { return i*(_n2+1)*(_n1+1) + j*(_n1+1) + k; }
    Uint t_idx(Uint i, Uint j, Uint k, Uint l) const { return i*_n2*_n1*6 + j*_n1*6 + k*6 + l; }

  protected:
    void buildBoundary(MultiGridCL* mgp) const;

  public:
    BBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&,
               Uint, Uint, Uint, Uint, Uint, Uint);

    virtual void
    build(MultiGridCL*) const;
};


class TetraBuilderCL : public MGBuilderCL
{
  private:
    const Ubyte rule_;

    const Point3DCL p0_;
    const Point3DCL p1_;
    const Point3DCL p2_;
    const Point3DCL p3_;

  protected:
    void buildBoundary(MultiGridCL* mgp) const;

  public:
    TetraBuilderCL(Ubyte rule);
    TetraBuilderCL(Ubyte rule, const Point3DCL& p0, const Point3DCL& p1,
                               const Point3DCL& p2, const Point3DCL& p3);

    static void
    BogoReMark(DROPS::MultiGridCL& mg, DROPS::Uint rule);
    virtual void
    build(MultiGridCL*) const;
};

/// \brief Create data structures for a tetrahedron-shaped multigrid
/** For detailed information see description of EmptyBrickBuilderCL*/
class EmptyTetraBuilderCL : public TetraBuilderCL
{
  private:
    Uint _numLevel;
    typedef TetraBuilderCL base_;

  public:
    EmptyTetraBuilderCL(Uint Level=1) : base_(0), _numLevel(Level) {}

    virtual void build( MultiGridCL*) const;
};

//--------------------------------------------------------------------
// Mesh-file-parser
//--------------------------------------------------------------------

// Reads a String as found in Mesh-Files:
// 1st eat whitespace and the following ".
// 2nd read string until second ". The second " is removed from the input, but
// not part of the string read.
// Usage: std::istream >> MeshStringCL >> std::string;
// TODO: Maybe handle "H\"allo" correctly.
class MeshStringCL
{
  private:
    std::istream* isp_;

  public:
    MeshStringCL()
        :isp_( 0) {}

    friend MeshStringCL&
    operator>>(std::istream&, MeshStringCL&);

    std::istream&
    operator>>(std::string&);

    static std::istream& // 2nd arg==true: Eat " ...... "
    SkipMeshString(std::istream&, bool= true); // 2nd arg==false: Eat .... "
};


// The headerof a mesh-file-section is typically a 5-tuple of hex-numbers, enclosed
// in parentheses.
typedef std::vector<Uint> HeaderInfoCL;

// Mesh files are organized in sections. See Appendix C of the TGrid User's Guide.
// Indices are 1-bsed.

// Node-Section
struct NodeSectionCL
{
    HeaderInfoCL headerinfo;
    std::vector<Point3DCL> point;
};

// All node sections of a file; only the section without node data is left out.
struct MeshNodeCL
{
    std::vector<NodeSectionCL> section;
    Uint num_expected;

    void
    Check( std::ostream* msg= &std::cout);
};

// Three node-ids, then id of right cell, then left cell. On boundaries, one of
// the neighbors may be zero, i. e. no cell.
typedef SArrayCL<Uint, 5> MFaceCL;

// Face-section
struct MFaceSectionCL
{
    HeaderInfoCL headerinfo;
    std::vector<MFaceCL> mface;
};

// All face sections of a file, apart from the section that only declares the
// total number of nodes.
struct MeshFaceCL
{
    std::vector<MFaceSectionCL> section;
    Uint num_expected;

    void // Check sanity of data and reoder sections.
    Check( std::ostream* msg= &std::cout);
    MFaceCL // 1-based index-access to mfaces
    operator[](Uint) const; // This could probably be sped up a lot if it should ever
                            // be neccessary. Call Check(), before you use this!!!
};

// Compares MFaceSectionCL based on their zone-id; used by algorithms in Check().
class FirstIndexLessCL : public std::binary_function<MFaceSectionCL, MFaceSectionCL, bool>
{
  public:
    bool operator () (const MFaceSectionCL& s0, const MFaceSectionCL& s1) const {
        return s0.headerinfo[1] < s1.headerinfo[1];
    }
};



// For one-based node-ids.
typedef SArrayCL<Uint, 4> CellCL;

// Cell-section
struct CellSectionCL
{
    HeaderInfoCL headerinfo;
//    std::vector<CellCL> cell;
};

// All cell-sections, apart from the one that only declares the total number of cells.
// Quite useless right now, as DROPS does not otherwise support volume-sections.
struct MeshCellCL
{
    std::vector<CellSectionCL> section;
    Uint num_expected;

    void
    Check( std::ostream* msg= &std::cout);
};

// Used to accumulate the faces and vertices of the cells, as this information
// is stored in the faces.
class HybridCellCL
{
  private:
    std::vector<FaceCL*> fp;
    std::vector<MFaceCL> mf;

  public:
    void
    push_back( FaceCL* fp_, MFaceCL mf_) {
        fp.push_back( fp_);
        mf.push_back( mf_);
    }
    std::vector<Uint> // Gather the vertices of this cell in consistent order from the faces.
    Vertices();
    std::vector<FaceCL*>
    Faces();
    FaceCL* // The face that has these nodes; compared as sets, not tuples.
    Face(Uint i, Uint j, Uint k);
    FaceCL* // Accessor to the faces; unordered
    face(Uint);
    void
    Check();
};

// Reads a mesh-file for FLUENT-UNS, FLUENT-RAMPANT and creates a multigrid on the
// call of build().
// The input-stream is only read, if build() is called.
// All temporary data is destroyed again on the exit of build() to conserve memory.
// We expect, that all simplexes are numbered consecutively starting with id 1.
// Currently, only a single node section and a single cell section are allowed. There
// may however be multiple face sections. Their handling could be extended easily to
// support multiple node and cell sections.
// Each boundary-face-section is translated into a Boundary-segment.
// TODO: It is utterly uncomfortable that build() must be const -- see all the mutables...
class ReadMeshBuilderCL : public MGBuilderCL
{
  private:
    static const char *SymbolicName_[];

    std::istream& f_;
    mutable std::vector<Uint> id_history_;
    mutable std::ostream* msg_;
    mutable std::vector<BndCondT> BC_;

    mutable std::map<Uint, BndIdxT> zone_id2bndidx_;


    bool // Find next (, that is not in a MeshString.
    NextSection() const;
         // Skip a section; for true: eat the leading ( of the section to be skipped,
    bool //                 for false: do not do this.
    SkipSection(bool= true) const;
         // Read a ( and the the following section-id, which is returned.
    Uint // The id is appended to the section-history.
    ReadId();
    HeaderInfoCL // Read ( hexnum0, hexnum1, .... )
    ReadHeaderInfoHex();
    void // Read a node-section
    ReadNode();
    void // Read a face-section
    ReadFace();
    void //Read a cell section; mostly useless.
    ReadCell();

    void // If not existent, add the boundary section given as the second argument.
    AddVertexBndDescription(VertexCL*, Uint) const;
         // Adds a boundary-idx to the edge given by the vertices. If the edge does not
    void // exist, it is created and recycled.
    CreateUpdateBndEdge(MultiGridCL::EdgeLevelCont&,
                        VertexCL*, VertexCL*, Uint) const;

    void // Deallocate memory of data-members
    Clear() const;

    static BndCondT MapBC( Uint gambit_bc); // map gambit bc to DROPS bc

  protected:
    mutable MeshNodeCL nodes_;
    mutable MeshFaceCL mfaces_;
    mutable MeshCellCL cells_;

    void buildBoundary(MultiGridCL*) const;
    // Actually creates the boundary-description into the MultiGrid; assumes
    // that the data file has been read.
    void buildBoundaryImp(MultiGridCL*) const;

    // Read a mesh-file; Uses above Read*-functions.
    void ReadFile();

  public:
    // Input stream, from which the mesh is read. Pass a pointer to an output stream,
    // e. g. msg= &std::cout, if you want to know, what happens during multigrid-construction.
    ReadMeshBuilderCL(std::istream& f, std::ostream* msg= 0);

    virtual void
    build(MultiGridCL*) const;

    static const char* // Symbolic section-names per TGrid User's Guide.
    Symbolic(Uint id);

    BndCondT GetBC( Uint i)                    const { return BC_[i]; }
    void     GetBC( std::vector<BndCondT>& BC) const { BC= BC_; }
};

/// \brief This class creates an empty domain out of a mesh-file
/** For detailed information see description of EmptyBrickBuilderCL*/
class EmptyReadMeshBuilderCL : public ReadMeshBuilderCL
{
  private:
    Uint numLevel_;
    typedef ReadMeshBuilderCL base_;

  public:
    EmptyReadMeshBuilderCL (std::istream& f, std::ostream* msg=0, Uint Level=1)
      : base_(f,msg) , numLevel_(Level) {}

    /// Just create the level views boundary-information
    virtual void build(MultiGridCL*) const;
};

/*******************************************************************
*   F I L E B U I L D E R  C L                                    *
*******************************************************************/

class FileBuilderCL : public MGBuilderCL
{
  private:
    // Path or File-Prefix
    std::string path_;

    MGBuilderCL* bndbuilder_;

    // Int <-> Add
    mutable std::map<size_t, VertexCL*> vertexAddressMap;
    mutable std::map<size_t, EdgeCL*>     edgeAddressMap;
    mutable std::map<size_t, FaceCL*>     faceAddressMap;
    mutable std::map<size_t, TetraCL*>   tetraAddressMap;

#ifdef _PAR
    /// \brief Read parallel information from a stream
    template <typename SimplexT>
    void ReadParInfo(std::istream&, SimplexT&) const;
#endif
    void BuildVerts   (MultiGridCL*) const;
    void BuildEdges   (MultiGridCL*) const;
    void BuildFacesI  (MultiGridCL*) const;
    void BuildTetras  (MultiGridCL*) const;
    void AddChildren  ()             const;
    void BuildFacesII (MultiGridCL*) const;
    void CheckFile( const std::ifstream& is) const;

  protected:
    void buildBoundary (MultiGridCL* mgp) const {bndbuilder_->buildBoundary(mgp);};

  public:
    FileBuilderCL (std::string path, MGBuilderCL* bndbuilder) : path_(path), bndbuilder_(bndbuilder) {};
    virtual void build(MultiGridCL*) const;
};


/*******************************************************************
*   M G S E R I A L I Z A T I O N   C L                           *
*******************************************************************/

class MGSerializationCL
{
  private:
    MultiGridCL& mg_;

    // Path or File-Prefix
    std::string  path_;

    // Addr <-> Int
    std::map<const EdgeCL*,   size_t>   edgeAddressMap;
    std::map<const VertexCL*, size_t> vertexAddressMap;
    std::map<const FaceCL*,   size_t>   faceAddressMap;
    std::map<const TetraCL*,  size_t>  tetraAddressMap;

    template<class itT, class T>
    void GetAddr (itT b, itT e, std::map<T, size_t> &m);

    void CreateAddrMaps ();

    // Writing-Routines
#ifdef _PAR
    /// \brief Write information about distribution
    template <typename SimplexT>
    void WriteParInfo(std::ofstream& os, const SimplexT&) const;
#endif
    void WriteEdges    ();
    void WriteFaces    ();
    void WriteVertices ();
    void WriteTetras   ();

    void CheckFile( const std::ofstream& os) const;

  public:
    MGSerializationCL (MultiGridCL& mg, std::string path) : mg_(mg), path_(path) {}
    void WriteMG ();
};

#ifdef _PAR

/// \brief Specialization of ReadParInfo for edges due to accMFR
template <>
void FileBuilderCL::ReadParInfo<EdgeCL>(std::istream&, EdgeCL&) const;

template <typename SimplexT>
  void FileBuilderCL::ReadParInfo(std::istream& is, SimplexT& s) const
/// read information from stream \a is, tell simplex these information. If
/// the simplex was a distributed object while serialization notify DDD about
/// the this distributed object
{
	PrioT prio;
    int numDist, distProc;
    GIDT oldGid;

    is >> prio >> numDist;
    s.SetPrio( prio);
    if ( numDist){
        is >> oldGid;
        for (int i=0; i<numDist; ++i){
            is >> distProc;
            Assert( distProc!=ProcCL::MyRank(), DROPSErrCL("FileBuilderCL::ReadParInfo: Cannot identify with myself"), DebugParallelC);
            DynamicDataInterfaceCL::IdentifyNumber( s.GetHdr(), distProc, oldGid);
        }
    }
}

template <typename SimplexT>
  void MGSerializationCL::WriteParInfo(std::ofstream& os, const SimplexT& s) const
{
    os << ' ' << s.GetPrio() << ' ' << s.GetNumDist()-1;
    if ( !s.IsLocal()){
        os << ' ' << s.GetGID();
        for (const int* proclist= s.GetProcList(); *proclist!=-1; proclist+=2)
            if ( *proclist!=ProcCL::MyRank())
                os << ' ' << (*proclist);
    }
}
#endif

} //end of namespace DROPS

#endif

//**************************************************************************
// File:    builder.h                                                      *
// Content: MGBuilderCL objects for some domains                           *
// Author:  Joerg Peters, Volker Reichelt, IGPM RWTH Aachen                *
// Version: 0.1                                                            *
// History: begin - October, 3 2000                                        *
//                                                                         *
// Remarks: We should use the const-qualifier to make it difficult to      *
//          accidentally change the multigrid structure from anywhere      *
//          outside of the multigrid algorithms.                           *
//          Thus the pointer to user data structures should probably be    *
//          a pointer to mutable.                                          *
//**************************************************************************


#ifndef DROPS_BUILDER_H
#define DROPS_BUILDER_H


#include "geom/multigrid.h"
#include <istream>

namespace DROPS
{

enum BndCondT
{
    DirBC= 0, Dir0BC= 2,         // (in)hom. Dirichlet boundary conditions
    NatBC= 1, Nat0BC= 3,         // (in)hom. natural   boundary conditions
    OutflowBC= 3, WallBC= 2,     // for convenience
    
    UndefinedBC_= -1             // ReadMeshBuilderCL: error, unknown bc
};


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

  public:
    BrickBuilderCL(const Point3DCL&, const Point3DCL&, const Point3DCL&, const Point3DCL&, Uint, Uint, Uint);

    virtual void
    build(MultiGridCL*) const;
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
    
  public:

    TetraBuilderCL(Ubyte rule);
    TetraBuilderCL(Ubyte rule, const Point3DCL& p0, const Point3DCL& p1,
                               const Point3DCL& p2, const Point3DCL& p3);

    static void
    BogoReMark(DROPS::MultiGridCL& mg, DROPS::Uint rule);
    virtual void
    build(MultiGridCL*) const;
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
    Check( std::ostream* msg= &std::cerr);
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
    Check( std::ostream* msg= &std::cerr);
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
// Quite useless right now, as DROPS does not otherwise support voume-sections.
struct MeshCellCL
{
    std::vector<CellSectionCL> section;
    Uint num_expected;

    void
    Check( std::ostream* msg= &std::cerr);
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

    mutable std::istream& f_;
    mutable std::vector<Uint> id_history_;
    mutable std::ostream* msg_;
    mutable std::vector<BndCondT> BC_;

    mutable MeshNodeCL nodes_;
    mutable MeshFaceCL mfaces_;
    mutable MeshCellCL cells_;

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
    void // Read a mesh-file; Uses above Read*-functions.
    ReadFile();

    void // If not existent, add the boundary section given as the second argument.
    AddVertexBndDescription(VertexCL*, Uint) const;
         // Adds a boundary-idx to the edge given by the vertices. If the edge does not
    void // exist, it is created and recycled.
    CreateUpdateBndEdge(MultiGridCL::EdgeLevelCont&,
                        VertexCL*, VertexCL*, Uint) const;

    void // Deallocate memory of data-members
    Clear() const;

    static BndCondT MapBC( Uint gambit_bc); // map gambit bc to DROPS bc

  public:
    // Input stream, from which the mesh is read. Pass a pointer to an output stream,
    // e. g. msg= &std::cerr, if you want to know, what happens during multigrid-construction.
    ReadMeshBuilderCL(std::istream& f, std::ostream* msg= 0);

    virtual void
    build(MultiGridCL*) const;

    static const char* // Symbolic section-names per TGrid User's Guide.
    Symbolic(Uint id);
    
    BndCondT GetBC( Uint i)                    const { return BC_[i]; }
    void     GetBC( std::vector<BndCondT>& BC) const { BC= BC_; }
};


} //end of namespace DROPS

#endif

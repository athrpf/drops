/// \file problem.h
/// \brief Classes that constitute a problem.
/// \author LNM RWTH Aachen: Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#ifndef DROPS_PROBLEM_H
#define DROPS_PROBLEM_H


#include "geom/multigrid.h"
#include "geom/builder.h"
#include "num/spmat.h"

#ifdef _PAR
#  include "parallel/parallel.h"   // for parallel reductions
#  include "parallel/interface.h"  // for accumulation of vectors
#endif

namespace DROPS
{

/// Prints a text-message describing the given boundary-condition.
void BndCondInfo (BndCondT, std::ostream&);

enum FiniteElementT
/// \brief enum for several FE types
///
/// Conventions:
/// - values < 128 are used for scalar FE
/// - values >= 128 are used for vector-valued FE,
///   the difference to the scalar FE counterpart should be 128
{
    P0_FE=0, P1_FE=1, P2_FE=2, P1Bubble_FE=3,  // for scalars
    P1D_FE=4, P1X_FE=5, P1IF_FE=6,
    vecP2_FE=130, vecP1Bubble_FE=131,          // for vectors
    UnknownFE_=-1
};

/// \brief For a given finite element type, this class describes how many degrees of freedom
/// are stored for each simplex-type on a single tetrahedron.
class FE_InfoCL
{
  protected:
    FiniteElementT fe_;    ///< FE type
    //@{
    /// \brief Number of unknowns on the simplex-type.
    Uint NumUnknownsVertex_;
    Uint NumUnknownsEdge_;
    Uint NumUnknownsFace_;
    Uint NumUnknownsTetra_;
    //@}

  public:
    FE_InfoCL( FiniteElementT fe) { SetFE(fe); }

    /// \brief Initialize with given FE-type \a fe
    void SetFE( FiniteElementT fe)
    {
        NumUnknownsVertex_= NumUnknownsEdge_= NumUnknownsFace_= NumUnknownsTetra_= 0;
        switch(fe_= fe) {
            case P0_FE:          NumUnknownsTetra_= 1; break;
            case P1_FE:
            case P1IF_FE:
            case P1X_FE:         NumUnknownsVertex_= 1; break;
            case P1Bubble_FE:    NumUnknownsVertex_= NumUnknownsTetra_= 1; break;
            case vecP1Bubble_FE: NumUnknownsVertex_= NumUnknownsTetra_= 3; break;
            case P1D_FE:         NumUnknownsFace_= 1; break;
            case P2_FE:          NumUnknownsVertex_= NumUnknownsEdge_= 1; break;
            case vecP2_FE:       NumUnknownsVertex_= NumUnknownsEdge_= 3; break;
            default:             throw DROPSErrCL("FE_InfoCL: unknown FE type");
        }
    }

    /// \brief Return enum code corresponding to FE-type
    FiniteElementT GetFE() const { return fe_; }
    /// \brief Returns true for XFEM
    bool IsExtended() const { return fe_==P1X_FE; }
    /// \brief Returns true for interface FE
    bool IsOnInterface() const { return fe_==P1IF_FE; }

    /// \brief Number of unknowns on the simplex-type
    //@{
    Uint NumUnknownsVertex() const { return NumUnknownsVertex_; }
    Uint NumUnknownsEdge()   const { return NumUnknownsEdge_; }
    Uint NumUnknownsFace()   const { return NumUnknownsFace_; }
    Uint NumUnknownsTetra()  const { return NumUnknownsTetra_; }
    template<class SimplexT> Uint GetNumUnknownsOnSimplex() const
    { throw DROPSErrCL("IdxDescCL::GetNumUnknowns: Unknown Simplex type"); }
    //@}
};

template<>
  inline Uint FE_InfoCL::GetNumUnknownsOnSimplex<VertexCL>() const { return NumUnknownsVertex(); }
template<>
  inline Uint FE_InfoCL::GetNumUnknownsOnSimplex<EdgeCL>()   const { return NumUnknownsEdge(); }
template<>
  inline Uint FE_InfoCL::GetNumUnknownsOnSimplex<FaceCL>()   const { return NumUnknownsFace(); }
template<>
  inline Uint FE_InfoCL::GetNumUnknownsOnSimplex<TetraCL>()  const { return NumUnknownsTetra(); }

class IdxDescCL;     // fwd decl
template<class T>
class VecDescBaseCL; // fwd decl

/// \brief The most widely used vector-description type; it uses a VectorCL
///     -object as Data.
typedef VecDescBaseCL<VectorCL> VecDescCL;


/// \brief Extended index for extended finite elements (XFEM).
///
/// Depending on the position of the zero-level of the levelset-function
/// additional DoFs are needed in the FE nodes. Internally, a std::vector
/// with a component for each FE node is kept. If the FE node has an
/// extended DoF, the index is stored in this component; otherwise,
/// NoIdx is stored.
class ExtIdxDescCL
{
    friend class IdxDescCL;
    typedef std::vector<IdxT> ExtendedIdxT;

  private:
    double       omit_bound_; ///< constant for stabilization of XFEM, controls omission of extended DoFs
    ExtendedIdxT Xidx_;       ///< vector entries store index of extended DoF (or NoIdx if FE node is not extended)
    ExtendedIdxT Xidx_old_;   ///< old extended index, used by member function Old2New(...)
#ifdef _PAR
    static IdxDescCL* current_Idx_;  ///< for DDD handlers
#endif

    ExtIdxDescCL( double omit_bound= 1./32. ) : omit_bound_( omit_bound ) {}

    /// \brief Update numbering of extended DoFs and return dimension of XFEM space (= # DoF + # extended DoF)
    ///
    /// Has to be called in two situations:
    /// - whenever level set function has changed to account for the moving interface (set \p NumberingChanged=false)
    /// - when numbering of index has changed, i.e. \p CreateNumbering was called before (set \p NumberingChanged=true)
    IdxT UpdateXNumbering( IdxDescCL*, const MultiGridCL&, const VecDescCL&, bool NumberingChanged= false );
    /// \brief Delete extended numbering
    void DeleteXNumbering() { Xidx_.resize(0); Xidx_old_.resize(0); }

  public:
    /// Get XFEM stabilization bound
    double GetBound() const { return omit_bound_; }
    /// Set XFEM stabilization bound
    void   SetBound( double omit_bound ) { omit_bound_= omit_bound; }
    /// Get extended index for DoF \p i
    IdxT  operator[]( const IdxT i ) const { return Xidx_[i]; }
    /// Set extended index for DoF \p i
    IdxT& operator[]( const IdxT i )       { return Xidx_[i]; }
    /// Returns number of DoFs for standard FE space
    IdxT  GetNumUnknownsStdFE()      const { return Xidx_.size(); }

    /// \brief Replace vector entries to account for different numbering of extended DoFs after call of UpdateXNumbering(...)
    void Old2New( VecDescCL* );
#ifdef _PAR
    /// \brief Gather xdof information on a vertex (for DDD)
    static int HandlerGatherUpdateXNumb ( OBJT objp, void* buf);
    /// \brief Scatter xdof information on a vertex (for DDD)
    static int HandlerScatterUpdateXNumb( OBJT objp, void* buf);
#endif
};

#ifdef _PAR
/// \name Wrapper for gathering and scattering data for ExtIdxDescCL::UpdateXNumbering
//@{
extern "C" inline int HandlerGatherUpdateXNumbC (OBJT objp, void* buf) { return ExtIdxDescCL::HandlerGatherUpdateXNumb ( objp, buf); }
extern "C" inline int HandlerScatterUpdateXNumbC(OBJT objp, void* buf) { return ExtIdxDescCL::HandlerScatterUpdateXNumb( objp, buf); }
//@}
class ExchangeCL;
#endif

/// \brief Mapping from the simplices in a triangulation to the components
///     of algebraic data-structures.
///
/// Internally, each object of type IdxDescCL has a unique index that is
/// used to access the unknown-indices that are stored in a helper class
/// (UnknownIdxCL and UnknownHandleCL) for each simplex. The unknown-indices
/// are allocated and numbered by using CreateNumbering.
class IdxDescCL: public FE_InfoCL
{
  private:
    static const Uint        InvalidIdx;   ///< Constant representing an invalid index.
    static std::vector<bool> IdxFree;      ///< Cache for unused indices; reduces memory-usage.

    Uint                     Idx_;         ///< The unique index.
    Uint                     TriangLevel_; ///< Triangulation of the index.
    IdxT                     NumUnknowns_; ///< total number of unknowns on the triangulation
    BndCondCL                Bnd_;         ///< boundary conditions
    match_fun                match_;       ///< matching function for periodic boundaries
    ExtIdxDescCL             extIdx_;      ///< extended index for XFEM
#ifdef _PAR
    ExchangeCL*              ex_;          ///< exchanging numerical data
#endif

    /// \brief Returns the lowest index that was not used and reserves it.
    Uint GetFreeIdx();
    /// \brief Number unknowns for standard FE.
    void CreateNumbStdFE( Uint level, MultiGridCL& mg);
    /// \brief Number unknowns on the vertices surrounding an interface.
    void CreateNumbOnInterface(Uint level, MultiGridCL& mg, const VecDescCL& ls, double omit_bound= -1./*default to using all dof*/);

  public:
    using FE_InfoCL::IsExtended;

    /// \brief The constructor uses the lowest available index for the
    ///     numbering. The triangulation level is initialized as 0.
    IdxDescCL( FiniteElementT fe= P1_FE, const BndCondCL& bnd= BndCondCL(0), match_fun match=0, double omit_bound=-99);
    /// \brief The copy will inherit the index number, whereas the index
    ///     of the original will be invalidated.
    IdxDescCL( const IdxDescCL& orig);
    /// \brief Frees the index, but does not invalidate the numbering on
    ///     the simplices. Call DeleteNumbering to do this.
    ~IdxDescCL();

    /// \brief Not implemented, as the private "Idx" should not be the
    ///     same for two different objects.
    IdxDescCL& operator= ( const IdxDescCL&);
    /// \brief Swaps the contents of obj and *this.
    void swap( IdxDescCL&);

    /// \brief Returns the number of the index. This can be used to access
    ///     the numbering on the simplices.
    Uint GetIdx() const {
        if (Idx_==InvalidIdx) throw DROPSErrCL("IdxDescCL::GetIdx: invalid index."
            " Probably using copy instead of original IdxDescCL-object.");
        return Idx_;
    }
    /// \brief Returns boundary condition
    BndCondCL GetBndInfo() const {return Bnd_;}
    /// \brief Returns extended index. Only makes sense for XFEM.
    const ExtIdxDescCL& GetXidx() const { return extIdx_; }
    /// \brief Returns extended index. Only makes sense for XFEM.
    ExtIdxDescCL&       GetXidx()       { return extIdx_; }
    /// \brief Triangulation of the index.
    Uint TriangLevel() const { return TriangLevel_; }
    /// \brief total number of unknowns on the triangulation
    IdxT NumUnknowns() const { return NumUnknowns_; }
    /// \brief Compare two IdxDescCL-objects. If a multigrid is given via mg, the
    ///     unknown-numbers on it are compared, too.
    static bool
    Equal(IdxDescCL& i, IdxDescCL& j, const MultiGridCL* mg= 0);

    /// \name Numbering
    /// \{
    /// \brief Used to number unknowns.
    void CreateNumbering( Uint level, MultiGridCL& mg, const VecDescCL* lsetp= 0);
    /// \brief Used to number unknowns and store boundary condition and matching function.
    void CreateNumbering( Uint level, MultiGridCL& mg, const BndCondCL& Bnd, match_fun match= 0, const VecDescCL* lsetp= 0)
    { Bnd_= Bnd; match_= match; CreateNumbering( level, mg, lsetp); }
    /// \brief Used to number unknowns taking boundary condition and matching function from \a baseIdx
    void CreateNumbering( Uint level, MultiGridCL& mg, const IdxDescCL& baseIdx, const VecDescCL* lsetp= 0)
    { Bnd_= baseIdx.Bnd_; match_= baseIdx.match_; CreateNumbering( level, mg, lsetp); }
    /// \brief Update numbering of extended DoFs.
    /// Has to be called whenever level set function has changed to account for the moving interface.
    void UpdateXNumbering( MultiGridCL& mg, const VecDescCL& lset);
    /// \brief Returns true, if XFEM is used and standard DoF \p dof is extended.
    bool IsExtended( IdxT dof) const
    { return IsExtended() ? extIdx_[dof] != NoIdx : false; }
    /// \brief Mark unknown-indices as invalid.
    void DeleteNumbering( MultiGridCL& mg);
    /// \}

#ifdef _PAR
    /// \brief Get a reference on the ExchangeCL
    ExchangeCL& GetEx() { return *ex_; }
    /// \brief Get a constant reference on the ExchangeCL
    const ExchangeCL& GetEx() const { return *ex_; }
    /// \brief get number of unknowns exclusively stored on this proc
    IdxT GetExclusiveNumUnknowns(const MultiGridCL&, int lvl=-1) const;
    /// \brief get global number of exclusive unknowns
    IdxT GetGlobalNumUnknowns(const MultiGridCL&, int lvl=-1) const;
#endif
};

/// \brief multilevel IdxDescCL
class MLIdxDescCL : public MLDataCL<IdxDescCL>
{
  public:
    MLIdxDescCL( FiniteElementT fe= P1_FE, size_t numLvl=1, const BndCondCL& bnd= BndCondCL(0), match_fun match=0, double omit_bound=1./32.)
    {
#ifdef _PAR
        if ( numLvl>1 )
            throw DROPSErrCL("MLIdxDescCL::MLIdxDescCL: No multilevel implemented in parDROPS, yet, sorry");
#endif
        for (size_t i=0; i< numLvl; ++i)
            this->push_back(IdxDescCL( fe, bnd, match, omit_bound));
    }

    void resize( size_t numLvl=1, FiniteElementT fe= P1_FE, const BndCondCL& bnd= BndCondCL(0), match_fun match=0, double omit_bound=1./32.)
    {
#ifdef _PAR
        if ( numLvl>1 )
            throw DROPSErrCL("MLIdxDescCL::resize: No multilevel implemented in parDROPS, yet, sorry");
#endif
        while (this->size() > numLvl)
            this->pop_back();
        while (this->size() < numLvl)
            this->push_back(IdxDescCL( fe, bnd, match, omit_bound));
    }

    /// \brief Returns the number of the index on the finest level.
    Uint GetIdx() const { return this->GetFinest().GetIdx();}

    /// \brief Triangulation of the index.
    Uint TriangLevel() const { return this->GetFinest().TriangLevel(); }

    /// \brief Returns true, if XFEM is used and standard DoF \p dof is extended.
    bool IsExtended( IdxT dof) const { return this->GetFinest().IsExtended( dof); }

    /// \brief total number of unknowns on the triangulation
    IdxT NumUnknowns() const { return this->GetFinest().NumUnknowns(); }

    /// \brief Number of unknowns on the simplex-type
    //@{
    Uint NumUnknownsVertex() const { return this->GetFinest().NumUnknownsVertex(); }
    Uint NumUnknownsEdge()   const { return this->GetFinest().NumUnknownsEdge(); }
    Uint NumUnknownsFace()   const { return this->GetFinest().NumUnknownsFace(); }
    Uint NumUnknownsTetra()  const { return this->GetFinest().NumUnknownsTetra(); }
    template<class SimplexT> Uint GetNumUnknownsOnSimplex() const
    { return this->GetFinest().GetNumUnknownsOnSimplex<SimplexT>(); }
    //@}

    /// \brief Initialize with given FE-type \a fe
    void SetFE( FiniteElementT fe )
    {
        for (MLIdxDescCL::iterator it = this->begin(); it != this->end(); ++it)
            it->SetFE(fe);
    }
    /// \name Numbering
    /// \{
    /// \brief Used to number unknowns on all levels.
    void CreateNumbering( size_t f_level, MultiGridCL& mg, const VecDescCL* lsetp= 0)
    {
        size_t lvl = ( this->size() <= f_level+1) ? f_level+1 - this->size(): 0;
        for (MLIdxDescCL::iterator it = this->begin(); it != this->end(); ++it)
        {
            it->CreateNumbering( lvl, mg, lsetp);
            if (lvl < f_level) ++lvl;
        }
    }
    /// \brief Used to number unknowns and store boundary condition and matching function on all levels.
    void CreateNumbering( size_t f_level, MultiGridCL& mg, const BndCondCL& Bnd, match_fun match= 0, const VecDescCL* lsetp= 0)
    {
        size_t lvl = ( this->size() <= f_level+1) ? f_level+1 - this->size(): 0;
        for (MLIdxDescCL::iterator it = this->begin(); it != this->end(); ++it)
        {
            it->CreateNumbering( lvl, mg, Bnd, match, lsetp);
            if (lvl < f_level) ++lvl;
        }
    }
    /// \brief Update numbering of extended DoFs on all levels.
    /// Has to be called whenever level set function has changed to account for the moving interface.
    void UpdateXNumbering( MultiGridCL& mg, const VecDescCL& lset)
    {
        for (MLIdxDescCL::iterator it = this->begin(); it != this->end(); ++it)
            it->UpdateXNumbering( mg, lset);
    }
    /// \brief Mark unknown-indices as invalid on all levels.
    void DeleteNumbering( MultiGridCL& mg)
    {
        for (MLIdxDescCL::iterator it = this->begin(); it != this->end(); ++it)
            it->DeleteNumbering( mg);
    }
    /// \}

#ifdef _PAR
    /// \brief Get a reference on the ExchangeCL (of the finest level)
    ExchangeCL& GetEx() { return this->GetFinest().GetEx(); }
    /// \brief Get a constant reference on the ExchangeCL (of the finest level)
    const ExchangeCL& GetEx() const { return this->GetFinest().GetEx(); }
    /// \brief get number of unknowns exclusively stored on this proc
    IdxT GetExclusiveNumUnknowns(const MultiGridCL& mg, int lvl=-1) const
    { return this->GetFinest().GetExclusiveNumUnknowns(mg, lvl); }
    /// \brief get global number of exclusive unknowns
    IdxT GetGlobalNumUnknowns(const MultiGridCL& mg, int lvl=-1) const
    { return this->GetFinest().GetGlobalNumUnknowns(mg, lvl); }
#endif
};

/// merges two p1-VectorCL into a p1x-VectorCL
void P1toP1X ( const IdxDescCL& xidx, VectorCL& p1x, const IdxDescCL& idx, const VectorCL& posPart, const VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg );

/// splits a p1x-VectorCL into two p1-VectorCL
void P1XtoP1 ( const IdxDescCL& xidx, const VectorCL& p1x, const IdxDescCL& idx, VectorCL& posPart, VectorCL& negPart, const VecDescCL& lset, const MultiGridCL& mg );

inline void
GetLocalNumbP1NoBnd(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx)
/// Copies P1-unknown-indices from idx on s into Numb; assumes that all
/// vertices have unknowns (NoBndDataCL and the like).
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i) {
        Numb[i]= s.GetVertex( i)->Unknowns( sys);
    }
}

inline void
GetLocalNumbP1DNoBnd(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx)
/// Copies P1D-unknown-indices from idx on s into Numb; assumes that all
/// faces have unknowns (NoBndDataCL and the like).
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i) {
        Numb[i]= s.GetFace( i)->Unknowns( sys);
    }
}

inline void
GetLocalNumbP2NoBnd(IdxT* Numb, const TetraCL& s, const IdxDescCL& idx)
/// Copies P2-unknown-indices from idx on s into Numb; assumes that all
/// vertices have unknowns (NoBndDataCL and the like).
{
    const Uint sys= idx.GetIdx();
    for (Uint i= 0; i < 4; ++i)
        Numb[i]= s.GetVertex( i)->Unknowns( sys);
    for(Uint i= 0; i < 6; ++i)
        Numb[i+4]= s.GetEdge( i)->Unknowns( sys);
}


/// \brief Collect indices of unknowns, boundary-segments and boundary
///     conditions on a tetrahedron.
///
/// This is convenient for discretisation of operators in the Setup-routines.
class LocalNumbP1CL
{
  public:
    /// \brief Field of unknown-indices; NoIdx, iff the degree of freedom lies
    /// on a boundary without unknowns. (Formerly called Numb.)
    IdxT     num   [4];
    /// \brief On boundaries, the number of the relevant BndSegDataCL-object
    /// in the corresponding BndDataCL-object, else NoBndC.
    BndIdxT  bndnum[4];
    /// \brief The relevant BndCondT, NoBC in the interior dofs.
    BndCondT bc    [4];

    /// \brief The default constructors leaves everything uninitialized.
    LocalNumbP1CL() {}
    /// \brief Read indices, boundary-segment numbers and boundary conditions
    ///     from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      LocalNumbP1CL(const TetraCL&, const IdxDescCL&, const BndDataT&);

    /// \brief Read indices, boundary-segment numbers and boundary conditions
    ///     from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      void
      assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd);

    /// \brief True, iff index i has a dof associated with it.
    bool WithUnknowns(IdxT i) const { return num[i] != NoIdx; }
};

/// \brief Collect indices of unknowns, boundary-segments and boundary
///     conditions on a tetrahedron.
///
/// This is convenient for discretisation of operators in the Setup-routines.
class LocalNumbP2CL
{
  public:
    /// \brief Field of unknown-indices; NoIdx, iff the degree of freedom lies
    /// on a boundary without unknowns. (Formerly called Numb.)
    IdxT     num   [10];
    /// \brief On boundaries, the number of the relevant BndSegDataCL-object
    /// in the corresponding BndDataCL-object, else NoBndC.
    BndIdxT  bndnum[10];
    /// \brief The relevant BndCondT, NoBC in the interior dofs.
    BndCondT bc    [10];

    /// \brief The default constructors leaves everything uninitialized.
    LocalNumbP2CL() {}
    /// \brief Read indices, boundary-segment numbers and boundary conditions
    ///     from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      LocalNumbP2CL(const TetraCL&, const IdxDescCL&, const BndDataT&);

    /// \brief Read indices, boundary-segment numbers and boundary conditions
    ///     from a tetrahedron and a BndDataCL-like object.
    template<class BndDataT>
      void
      assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd);

    /// \brief True, iff index i has a dof associated with it.
    bool WithUnknowns(IdxT i) const { return num[i] != NoIdx; }
};


/// \brief A numerical vector together with an IdxDescCL -object,
///     that couples it to simplices in a multigrid.
template<class T>
class VecDescBaseCL
{
  public:
    /// \brief The type of the numerical vector.
    typedef T DataType;

    /// \brief The default-constructor creates an empty vector and sets RowIdx to 0.
    VecDescBaseCL()
        :RowIdx(0) {}
    /// \brief Initialize RowIdx with idx and contruct Data with the given size.
    VecDescBaseCL( IdxDescCL* idx) { SetIdx( idx); }
    VecDescBaseCL( MLIdxDescCL* idx) { SetIdx( &(idx->GetFinest()) ); }


    IdxDescCL* RowIdx; ///< Pointer to the index-description used for Data.
    DataType  Data;    ///< The numerical data.

    /// \brief The triangulation-level of the index.
    Uint GetLevel() const { return RowIdx->TriangLevel(); }
    /// \brief Use a new index for accessing the components.
    void SetIdx(IdxDescCL*);
    void SetIdx(MLIdxDescCL* idx) {SetIdx( &( idx->GetFinest()) );}
    /// \brief Resize a vector according to RowIdx.
    void Clear();
    /// \brief Empty Data and set RowIdx to 0.
    void Reset();

    /// \brief Write Data on a stream
    void Write(std::ostream&, bool binary=false) const;
    /// \brief Read Data from a stream
    void Read(std::istream&, bool binary=false);
};


/// \brief A sparse matrix together with two IdxDescCL -objects,
///     that couple the row- and column- indices to simplices in a
///     multigrid.
template<typename MatT, typename IdxDescT>
class MatDescBaseCL
{
  public:
    /// \brief The type of the matrix.
    typedef MatT DataType;

    /// \brief The default-constructor creates an empty matrix and sets
    /// the index-pointers to 0.
    MatDescBaseCL()
        : RowIdx( 0), ColIdx( 0) {}
    /// \brief Initialize RowIdx an ColIdx; Data is still default-constructed.
    MatDescBaseCL( IdxDescT* r, IdxDescT* c) { SetIdx( r, c); }

    IdxDescT* RowIdx; ///< Pointer to the index-description used for row-indices.
    IdxDescT* ColIdx; ///< Pointer to the index-description used for column-indices.
    DataType  Data;   ///< The numerical data.

    /// \brief The triangulation-level of the row-index.
    Uint GetRowLevel() const { return RowIdx->TriangLevel(); }
    /// \brief The triangulation-level of the column-index.
    Uint GetColLevel() const { return ColIdx->TriangLevel(); }

    /// \brief Use a new index for accessing the components.
    void SetIdx( IdxDescT*, IdxDescT*);
    /// \brief Empty Data and set the index-pointers to 0.
    void Reset();
};

typedef MatDescBaseCL<MatrixCL,IdxDescCL>     MatDescCL;
typedef MatDescBaseCL<MLMatrixCL,MLIdxDescCL> MLMatDescCL;



/// \brief This class contains the main constituents of a forward problem
///     with a PDE.
///
/// \param Coeff   coefficients of the underlying PDE
/// \param BndData boundary conditions
/// \param ExCL    class to manage exchange of numerical values over processors
///
/// \todo Probably we should not copy CoeffCL and BndDataCL.
template <class Coeff, class BndData>
class ProblemCL
{
  public:
    typedef Coeff    CoeffCL;
    typedef BndData  BndDataCL;

  protected:
    bool         _myMG;
    MultiGridCL& _MG;         ///< The multigrid.
    CoeffCL      _Coeff;      ///< Right-hand-side, coefficients of the PDE.
    BndDataCL    _BndData;    ///< boundary-conditions

  public:
    /// \brief The multigrid constructed from mgbuilder will be destroyed if this variable leaves its scope.
    ProblemCL(const MGBuilderCL& mgbuilder, const CoeffCL& coeff, const BndDataCL& bnddata)
        : _myMG( true), _MG( *new MultiGridCL( mgbuilder)), _Coeff( coeff), _BndData( bnddata) {}
    /// \brief The multigrid mg will be left alone if this variable leaves its scope.
    ProblemCL(MultiGridCL& mg, const CoeffCL& coeff, const BndDataCL& bnddata)
        : _myMG( false), _MG( mg), _Coeff( coeff), _BndData( bnddata) {}
    ~ProblemCL() { if (_myMG) delete &_MG; }

    MultiGridCL&       GetMG()            { return _MG; }
    const MultiGridCL& GetMG()      const { return _MG; }
    const CoeffCL&     GetCoeff()   const { return _Coeff; }
    const BndDataCL&   GetBndData() const { return _BndData; }
};


template<class BndDataT>
  LocalNumbP1CL::LocalNumbP1CL (const TetraCL& s, const IdxDescCL& idx,
      const BndDataT& bnd)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL  -object to be used.
/// \param bnd The BndDataCL -like-object, from which boundary-segment-numbers are used.
{
    this->assign( s, idx, bnd);
}

template<class BndDataT>
  void
  LocalNumbP1CL::assign (const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL -object to be used.
/// \param bnd The BndDataCL -like object, from which boundary-segment-numbers are used.
{
    BndIdxT bidx= 0;
    const Uint sys= idx.GetIdx();

    for (Uint i= 0; i < NumVertsC; ++i)
        if (NoBC == (bc[i]= bnd.GetBC( *s.GetVertex( i), bidx))) {
            bndnum[i]= NoBndC;
            num[i]= s.GetVertex( i)->Unknowns( sys);
        }
        else {
            bndnum[i]= bidx;
            num[i]= (bnd.GetBndSeg( bidx).WithUnknowns())
                ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
        }
}

template<class BndDataT>
  LocalNumbP2CL::LocalNumbP2CL(const TetraCL& s, const IdxDescCL& idx,
      const BndDataT& bnd)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL  -object to be used.
/// \param bnd The BndDataCL -like-object, from which boundary-segment-numbers are used.
{
    this->assign( s, idx, bnd);
}

template<class BndDataT>
  void
  LocalNumbP2CL::assign(const TetraCL& s, const IdxDescCL& idx, const BndDataT& bnd)
/// \param s The tet, from which index-numbers are read.
/// \param idx The IdxDescCL -object to be used.
/// \param bnd The BndDataCL -like object, from which boundary-segment-numbers are used.
{
    BndIdxT bidx= 0;
    const Uint sys= idx.GetIdx();

    for (Uint i= 0; i < NumVertsC; ++i)
        if (NoBC == (bc[i]= bnd.GetBC( *s.GetVertex( i), bidx))) {
            bndnum[i]= NoBndC;
            num[i]= s.GetVertex( i)->Unknowns( sys);
        }
        else {
            bndnum[i]= bidx;
            num[i]= (bnd.GetBndSeg( bidx).WithUnknowns())
                ? s.GetVertex( i)->Unknowns( sys) : NoIdx;
        }
    for (Uint i= 0; i< NumEdgesC; ++i)
        if (NoBC == (bc[i+NumVertsC]= bnd.GetBC( *s.GetEdge( i), bidx))) {
            bndnum[i+NumVertsC]= NoBndC;
            num[i+NumVertsC]= s.GetEdge( i)->Unknowns( sys);
        }
        else {
            bndnum[i+NumVertsC]= bidx;
            num[i+NumVertsC]= (bnd.GetBndSeg( bidx).WithUnknowns())
                ? s.GetEdge( i)->Unknowns( sys) : NoIdx;
        }
}

template<class T>
void VecDescBaseCL<T>::SetIdx(IdxDescCL* idx)
/// Prepares the vector for usage with a new index-object for
/// its components. The vector is resized to size 0 and
/// then to the new size.
{
    RowIdx = idx;
    Data.resize(0);
    Data.resize(idx->NumUnknowns());
}

template<class T>
void VecDescBaseCL<T>::Clear()
/// The vector is resized to size 0 and then resized to the size given
/// by RowIdx.
{
    Data.resize(0);
    Data.resize(RowIdx->NumUnknowns());
}

template<class T>
void VecDescBaseCL<T>::Reset()
/// Sets RowIdx to 0 and resizes the vector to size 0.
{
    RowIdx = 0;
    Data.resize(0);
}

template<class T>
void VecDescBaseCL<T>::Write(std::ostream& os, bool binary) const
/// Writes numerical data on a stream, which can be read by VecDescBaseCL::Read
/// \param os where to put the data
/// \param binary write out data in binary format
{
    if (binary){
        size_t numUnk= Data.size();
        os.write( (char*)(&(numUnk)), sizeof(size_t));
        os.write( (char*)Addr(Data), sizeof(typename T::value_type)*Data.size());
    }

    else {
        // Write out all digits
        os.precision(std::numeric_limits<typename T::value_type>::digits10);
        os << Data;
    }
}

template<class T>
void VecDescBaseCL<T>::Read(std::istream& is, bool binary)
/// Read data from stream \a is, which should have been created
/// by VecDescBaseCL::Write.
/// \param is where to read the data
/// \param binary read data in binary format
/// \pre CreateNumbering for RowIdx and SetIdx must have been
///      called
{
    // read number of unknowns for error checking
    size_t readUnk= 0;
    if (binary)
        is.read( (char*)(&readUnk), sizeof(size_t));
    else{
        is >> readUnk;
        is.seekg(0);    // rewind
    }

    // Check if number of unknowns is correct
    if (   ( !RowIdx->IsExtended() && RowIdx->NumUnknowns()!=readUnk)
         ||( RowIdx->IsExtended()  && RowIdx->GetXidx().GetNumUnknownsStdFE()!=readUnk ) )
    {
        throw DROPSErrCL("VecDescBaseCL::Read: Number of Unknowns does not match, wrong FE-type?");
    }

    // read data
    if (binary)
        is.read( (char*)(Addr(Data)), sizeof(typename T::value_type)*readUnk);
    else
        in(is, Data);
}

template<typename MatT, typename IdxDescT>
void MatDescBaseCL<MatT, IdxDescT>::SetIdx( IdxDescT* row, IdxDescT* col)
/// Prepares the matrix for usage with new index-objects for
/// its components. As constructing sparse matrices is fairly involved,
/// this routine does not modify Data. SparseMatBuilderCL should be used
/// to do this.
{
    RowIdx= row;
    ColIdx= col;
}

template<typename MatT, typename IdxDescT>
void MatDescBaseCL<MatT, IdxDescT>::Reset()
/// Sets the index-pointers to 0 and clears the matrix.
{
    RowIdx = 0;
    ColIdx = 0;
    Data.clear();
}

/// \name Helper routines to number unknowns.
/// These functions should not be used directly. CreateNumb is much more
/// comfortable and as efficient.
///
/// These functions allocate memory for the Unknown-indices in system
/// idx on all simplices of the indicated type between begin and end.
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is \# Unknowns+stride.
/// Simplices on Dirichlet boundaries are skipped.
/// \{
template<class SimplexT>
void CreateNumbOnSimplex( const Uint idx, IdxT& counter, Uint stride,
                         const ptr_iter<SimplexT>& begin,
                         const ptr_iter<SimplexT>& end,
                         const BndCondCL& Bnd)
{
    if (stride == 0) return;
    for (ptr_iter<SimplexT> it= begin; it != end; ++it)
    {
        if ( !Bnd.IsOnDirBnd( *it) )
        {
#ifdef _PAR
            if (it->IsMarkedForRemovement())
                continue;
#endif
            it->Unknowns.Prepare( idx);
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
}

void CreateNumbOnTetra( const Uint idx, IdxT& counter, Uint stride,
                        const MultiGridCL::TriangTetraIteratorCL& begin,
                        const MultiGridCL::TriangTetraIteratorCL& end);
/// \}


/// \brief Mark unknown-indices as invalid.
///
/// This routine writes NoIdx as unknown-index for all indices of the
/// given system.
/// \note Currently, there is no way to compactify the memory used by the
///     UnknownHandleCL -objects in the simplices, as we never had a high
///     variation in the number of allocated indices. Adding such a pass
///     is probably not hard and defered until the need for one arises.
/// \param idx The system-number, as returned by IdxDescCL::GetIdx(),
///     to be invalidated.
/// \param begin The beginning of the sequence of simplices, on which the
///     system has indices.
/// \param end The end of the sequence of simplices, on which the
///     system has indices.
///
/// \note To be sure that all indices are invalidated, one can use the
///     iterators returned by, e.g., MultiGridCL::GetAllVertexBegin,
///     MultiGridCL::GetAllVertexEnd.
template <class Iter>
inline void
DeleteNumbOnSimplex( Uint idx, const Iter& begin, const Iter& end)
{

    for (Iter it=begin; it!=end; ++it)
        if (it->Unknowns.Exist() && it->Unknowns.Exist( idx) )
            it->Unknowns.Invalidate( idx);
}


/// \name Helper routine to number unknowns on vertices, edges, faces on periodic boundaries.
/// This function should not be used directly. CreateNumb is much more
/// comfortable and as efficient.
///
/// For interior simplices and boundary-conditions other than Per1BC and Per2BC,
/// nothing unusual happens.
/// The matching works as follows:
///     - Simplices with Per1BC are memoized in list l1, those with Per2BC in list l2.
///     - The unknowns in l1 are numbered.
///     - Each element of l2 is matched in l1 via the matching function and inherits the
///       indices from its l1-counterpart.
/// \{
template<class SimplexT>
void CreatePeriodicNumbOnSimplex( const Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const ptr_iter<SimplexT>& begin,
                        const ptr_iter<SimplexT>& end,
                        const BndCondCL& Bnd)
{
    if (stride == 0) return;

    typedef std::list<SimplexT*> psetT;
    typedef typename std::list<SimplexT*>::iterator psetIterT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (ptr_iter<SimplexT> it= begin; it!=end; ++it)
    {
        if ( Bnd.IsOnDirBnd( *it) ) continue;
        it->Unknowns.Prepare( idx);
        if (Bnd.IsOnPerBnd( *it))
        {
            if (Bnd.GetBC( *it)==Per1BC)
            {
                s1.push_back( &*it);
                it->Unknowns( idx)= counter;
                counter+= stride;
            }
            else
                s2.push_back( &*it);
        }
        else
        {
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
    // now we have s1.size() <= s2.size()
    // match objects in s1 and s2
    for (psetIterT it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetIterT it2= s2.begin(), end2= s2.end(); it2!=end2;)
            if (match( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2
                s2.erase( it2++);
            }
            else it2++;
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnSimplex: Periodic boundaries do not match!");
}
/// \}

} // end of namespace DROPS


#endif

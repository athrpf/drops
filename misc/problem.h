/// \file
/// \brief Classes that constitute a problem.

#ifndef DROPS_PROBLEM_H
#define DROPS_PROBLEM_H


#include "geom/multigrid.h"
#include "geom/builder.h"
#include "num/spmat.h"

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
    P1D_FE=4, P1X_FE=5,
                      vecP2_FE= 130            // for vectors
};


/// \brief Mapping from the simplices in a triangulation to the components
///     of algebraic data-structures.
///
/// This class describes how many unknowns are reserved for each
/// simplex in a given triangulation. The number of unknowns can be
/// given separately for each simplex-type.  Note, that no memory for the
/// numbers or unknowns is allocated. This must be done in another step,
/// e.g. by calling CreateNumb.

/// Internally, each object of type IdxDescCL has a unique index that is
/// used to access the unknown-indices that are stored in a helper class
/// (UnknownIdxCL and UnknownHandleCL) for each simplex.
class IdxDescCL
{
  private:
    static const Uint        InvalidIdx; ///< Constant representing an invalid index.
    static std::vector<bool> IdxFree;    ///< Cache for unused indices; reduces memory-usage.
    Uint                     Idx;        ///< The unique index.

    /// \brief Returns the lowest index that was not used and reserves it.
    Uint GetFreeIdx();

  public:
    Uint TriangLevel;        ///< Triangulation of the index.

    //@{
    /// \brief Number of unknowns on the simplex-type.
    Uint NumUnknownsVertex;
    Uint NumUnknownsEdge;
    Uint NumUnknownsFace;
    Uint NumUnknownsTetra;
    IdxT NumUnknowns;
    //@}

    /// \brief The constructor uses the lowest available index for the
    ///     numbering. The triangulation level must be set separately.
    IdxDescCL( Uint unkVertex= 0, Uint unkEdge= 0, Uint unkFace= 0, Uint unkTetra= 0)
      : Idx( GetFreeIdx()), NumUnknownsVertex( unkVertex), NumUnknownsEdge( unkEdge),
        NumUnknownsFace( unkFace), NumUnknownsTetra( unkTetra), NumUnknowns( 0) {}
    explicit IdxDescCL( FiniteElementT fe)
      : Idx( GetFreeIdx()), NumUnknownsVertex( 0), NumUnknownsEdge( 0),
        NumUnknownsFace( 0), NumUnknownsTetra( 0), NumUnknowns( 0)
    {
        switch(fe) {
            case P0_FE:       NumUnknownsTetra= 1; break;
            case P1_FE:
            case P1X_FE:      NumUnknownsVertex= 1; break;
            case P1Bubble_FE: NumUnknownsVertex= NumUnknownsTetra= 1; break;
            case P1D_FE:      NumUnknownsFace= 1; break;
            case P2_FE:       NumUnknownsVertex= NumUnknownsEdge= 1; break;
            case vecP2_FE:    NumUnknownsVertex= NumUnknownsEdge= 3; break;
        }
    }
    /// \brief The copy will inherit the index number, whereas the index
    ///     of the original will be invalidated.
    IdxDescCL( const IdxDescCL& orig);
    /// \brief Frees the index, but does not invalidate the numbering on
    ///     the simplices. Call DeleteNumbering to do this.
    ~IdxDescCL() { if (Idx!=InvalidIdx) IdxFree[Idx]= true; }

    /// \brief Not implemented, as the private "Idx" should not be the
    ///     same for two different objects.
    IdxDescCL& operator= ( const IdxDescCL&);
    /// \brief Swaps the contents of obj and *this.
    void swap( IdxDescCL&);

    /// \brief Set the number of unknowns per simplex-type for this index.
    void Set( Uint unkVertex, Uint unkEdge= 0, Uint unkFace= 0, Uint unkTetra= 0);
    /// \brief Returns the number of the index. This can be used to access
    ///     the numbering on the simplices.
    Uint GetIdx() const {
        if (Idx==InvalidIdx) throw DROPSErrCL("IdxDescCL::GetIdx: invalid index."
            " Probably using copy instead of original IdxDescCL-object.");
        return Idx;
    }

    /// \brief Compare two IdxDescCL-objects. If a multigrid is given via mg, the
    ///     unknown-numbers on it are compared, too.
    static bool
    Equal(IdxDescCL& i, IdxDescCL& j, const MultiGridCL* mg= 0);
};


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

    IdxDescCL* RowIdx; ///< Pointer to the index-description used for Data.
    DataType  Data;    ///< The numerical data.

    /// \brief The triangulation-level of the index.
    Uint GetLevel() const { return RowIdx->TriangLevel; }
    /// \brief Use a new index for accessing the components.
    void SetIdx(IdxDescCL*);
    /// \brief Resize a vector according to RowIdx.
    void Clear();
    /// \brief Empty Data and set RowIdx to 0.
    void Reset();
};


/// \brief The most widely used vector-description type; it uses a VectorCL
///     -object as Data.
typedef VecDescBaseCL<VectorCL> VecDescCL;


/// \brief A sparse matrix together with two IdxDescCL -objects,
///     that couple the row- and column- indices to simplices in a
///     multigrid.
class MatDescCL
{
  public:
    /// \brief The type of the matrix.
    typedef MatrixCL DataType;

    /// \brief The default-constructor creates an empty matrix and sets
    /// the index-pointers to 0.
    MatDescCL()
        :RowIdx(0), ColIdx(0) {}
    /// \brief Initialize RowIdx an ColIdx; Data is still default-constructed.
    MatDescCL(IdxDescCL* r, IdxDescCL* c) { SetIdx( r, c); }

    IdxDescCL* RowIdx; ///< Pointer to the index-description used for row-indices.
    IdxDescCL* ColIdx; ///< Pointer to the index-description used for column-indices.
    DataType  Data; ///< The numerical data.

    /// \brief The triangulation-level of the row-index.
    Uint GetRowLevel() const { return RowIdx->TriangLevel; }
    /// \brief The triangulation-level of the column-index.
    Uint GetColLevel() const { return ColIdx->TriangLevel; }

    /// \brief Use a new index for accessing the components.
    void SetIdx(IdxDescCL*, IdxDescCL*);
    /// \brief Empty Data and set the index-pointers to 0.
    void Reset();
};


/// \brief This class contains the main constituents of a forward problem
///     with a PDE.
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
    Data.resize(idx->NumUnknowns);
}

template<class T>
void VecDescBaseCL<T>::Clear()
/// The vector is resized to size 0 and then resized to the size given
/// by RowIdx.
{
    Data.resize(0);
    Data.resize(RowIdx->NumUnknowns);
}

template<class T>
void VecDescBaseCL<T>::Reset()
/// Sets RowIdx to 0 and resizes the vector to size 0.
{
    RowIdx = 0;
    Data.resize(0);
}


/// \name Routines to number unknowns.
/// These functions should not be used directly. CreateNumb is much more
/// comfortable and as efficient.
///
/// These functions allocate memory for the Unknown-indices in system
/// idx on all simplices of the indicated type between begin and end.
/// The first number used is the initial value of counter, the next
/// numbers are counter+stride, counter+2*stride, and so on.
/// Upon return, counter contains the first number, that was not used,
/// that is #Unknowns+stride.
/// Simplices on Dirichlet boundaries are skipped.
/// \{
template<class BndDataT>
void CreateNumbOnVertex( const Uint idx, IdxT& counter, Uint stride,
                         const MultiGridCL::TriangVertexIteratorCL& begin,
                         const MultiGridCL::TriangVertexIteratorCL& end,
                         const BndDataT& Bnd)
{
    if (stride == 0) return;
    for (MultiGridCL::TriangVertexIteratorCL it= begin; it != end; ++it)
    {
        if ( !Bnd.IsOnDirBnd( *it) )
        {
            it->Unknowns.Prepare( idx);
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
}

template<class BndDataT>
void CreateNumbOnEdge( const Uint idx, IdxT& counter, Uint stride,
                       const MultiGridCL::TriangEdgeIteratorCL& begin,
                       const MultiGridCL::TriangEdgeIteratorCL& end,
                       const BndDataT& Bnd)
{
    if (stride == 0) return;
    for (MultiGridCL::TriangEdgeIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {
            it->Unknowns.Prepare( idx);
            it->Unknowns( idx)= counter;
            counter+= stride;
        }
    }
}

template<class BndDataT>
void CreateNumbOnFace( const Uint idx, IdxT& counter, Uint stride,
                       const MultiGridCL::TriangFaceIteratorCL& begin,
                       const MultiGridCL::TriangFaceIteratorCL& end,
                       const BndDataT& Bnd)
{
    if (stride == 0) return;
    for (MultiGridCL::TriangFaceIteratorCL it=begin; it!=end; ++it)
    {
        if ( !Bnd.IsOnDirBnd(*it) )
        {
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


/// \brief Type of functions used to identify points on periodic boundaries,
///     that share the same dof.
typedef bool (*match_fun)(const Point3DCL&, const Point3DCL&);


/// \name Routines to number unknowns on periodic boundaries.
/// These functions should not be used directly. CreateNumb is much more
/// comfortable and as efficient.
///
/// For interior simplices and boundary-conditions other than Per1BC and Per2BC,
/// nothing unusual happens.
/// The matching works as follows:
///     - Simplices with Per1BC are memoized in list l1, those with Per2BC in list l2.
///     - The unknowns in l1 are numbered.
///     - Each element of l2 is matched in l1 via the matching function and inherits the
///       indices from its l1-counterpart.
///     .
/// \todo In the presence of adaptive refinement the multigrid may be
///     such, that there is
///     no correspondence between the simplices on boundary-parts that are
///     identified for the boundary condition. However, enforcing such
///     a symmetry requires a modification of the refinement algorithm
///     which is probably not a trivial exercise.
/// \{
template<class BndDataT>
void CreatePeriodicNumbOnVertex( Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const MultiGridCL::TriangVertexIteratorCL& begin,
                        const MultiGridCL::TriangVertexIteratorCL& end,
                        const BndDataT& Bnd)
{
    if (stride == 0) return;

    typedef std::list<VertexCL*> psetT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (MultiGridCL::TriangVertexIteratorCL it= begin; it!=end; ++it)
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
    if (s1.size() != s2.size())
        throw DROPSErrCL( "CreatePeriodicNumbOnVertex: Periodic boundaries do not match!");
    // match objects in s1 and s2
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; ++it2)
            if (match( (*it1)->GetCoord(), (*it2)->GetCoord()))
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2 and stop search
                s2.erase( it2);
                break;
            }
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnVertex: Periodic boundaries do not match!");
}

template<class BndDataT>
void CreatePeriodicNumbOnEdge( Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const MultiGridCL::TriangEdgeIteratorCL& begin,
                        const MultiGridCL::TriangEdgeIteratorCL& end,
                        const BndDataT& Bnd)
{
    if (stride == 0) return;

    typedef std::list<EdgeCL*> psetT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (MultiGridCL::TriangEdgeIteratorCL it= begin; it!=end; ++it)
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
    if (s1.size() != s2.size())
        throw DROPSErrCL( "CreatePeriodicNumbOnEdge: Periodic boundaries do not match!");
    // match objects in s1 and s2
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; ++it2)
            if (match( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2 and stop search
                s2.erase( it2);
                break;
            }
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnEdge: Periodic boundaries do not match!");
}

template<class BndDataT>
void CreatePeriodicNumbOnFace( Uint idx, IdxT& counter, Uint stride, match_fun match,
                        const MultiGridCL::TriangFaceIteratorCL& begin,
                        const MultiGridCL::TriangFaceIteratorCL& end,
                        const BndDataT& Bnd)
{
    if (stride == 0) return;

    typedef std::list<FaceCL*> psetT;
    psetT s1, s2;
    // create numbering for all objects (skipping Dir bnds) except those on Per2 bnds.
    // collect all objects on Per1/Per2 bnds in s1, s2 resp.
    for (MultiGridCL::TriangFaceIteratorCL it= begin; it!=end; ++it)
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
    if (s1.size() != s2.size())
        throw DROPSErrCL( "CreatePeriodicNumbOnFace: Periodic boundaries do not match!");
    // match objects in s1 and s2
    for (psetT::iterator it1= s1.begin(), end1= s1.end(); it1!=end1; ++it1)
    {
        // search corresponding object in s2
        for (psetT::iterator it2= s2.begin(), end2= s2.end(); it2!=end2; ++it2)
            if (match( GetBaryCenter( **it1), GetBaryCenter( **it2)) )
            {
                // it2 gets same number as it1
                (*it2)->Unknowns( idx)= (*it1)->Unknowns( idx);
                // remove it2 from s2 and stop search
                s2.erase( it2);
                break;
            }
    }
    if (!s2.empty())
        throw DROPSErrCL( "CreatePeriodicNumbOnFace: Periodic boundaries do not match!");
}
/// \}


/// \brief Used to number unknowns depending on the index idx.
///
/// If a matching function is given, numbering on periodic boundaries
/// is performed, too.
/// Memory for the Unknown-Indices on TriangLevel level is allocated
/// and the unknowns are numbered.
/// \param level Level of the triangulation to use.
/// \param idx The index-description to be used.
/// \param mg The multigrid to be operated on.
/// \param Bnd A BndDataCL -like object used to determine the boundary-condition of boundary simplices.
/// \param match An optional matching function. It is only neccessary, if periodic boundaries are used.
/// \pre idx.NumUnknownsVertex etc. must be set.
/// \post The function sets idx.TriangLevel to level and sets idx.NumUnknowns.
template<class BndDataT>
void CreateNumb(Uint level, IdxDescCL& idx, MultiGridCL& mg, const BndDataT& Bnd, match_fun match= 0)
{
    // set up the index description
    idx.TriangLevel = level;
    idx.NumUnknowns = 0;

    const Uint idxnum= idx.GetIdx();

    // allocate space for indices; number unknowns in TriangLevel level
    if (match)
    {
        if (idx.NumUnknownsVertex)
            CreatePeriodicNumbOnVertex( idxnum, idx.NumUnknowns, idx.NumUnknownsVertex, match,
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd);
        if (idx.NumUnknownsEdge)
            CreatePeriodicNumbOnEdge( idxnum, idx.NumUnknowns, idx.NumUnknownsEdge, match,
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd);
        if (idx.NumUnknownsFace)
            CreatePeriodicNumbOnFace( idxnum, idx.NumUnknowns, idx.NumUnknownsFace, match,
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd);
        if (idx.NumUnknownsTetra)
            CreateNumbOnTetra( idxnum, idx.NumUnknowns, idx.NumUnknownsTetra,
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level));
    }
    else
    {
        if (idx.NumUnknownsVertex)
            CreateNumbOnVertex( idxnum, idx.NumUnknowns, idx.NumUnknownsVertex,
                mg.GetTriangVertexBegin(level), mg.GetTriangVertexEnd(level), Bnd);
        if (idx.NumUnknownsEdge)
            CreateNumbOnEdge( idxnum, idx.NumUnknowns, idx.NumUnknownsEdge,
                mg.GetTriangEdgeBegin(level), mg.GetTriangEdgeEnd(level), Bnd);
        if (idx.NumUnknownsFace)
            CreateNumbOnFace( idxnum, idx.NumUnknowns, idx.NumUnknownsFace,
                mg.GetTriangFaceBegin(level), mg.GetTriangFaceEnd(level), Bnd);
        if (idx.NumUnknownsTetra)
            CreateNumbOnTetra( idxnum, idx.NumUnknowns, idx.NumUnknownsTetra,
                mg.GetTriangTetraBegin(level), mg.GetTriangTetraEnd(level));
    }
}

/// \brief Mark unknown-indices as invalid for given index-description.
///
/// This routine writes NoIdx as unknown-index for all indices of the
/// given index-description. idx.NumUnknowns will be set to zero.
/// \param idx The index-description to be used.
/// \param MG The multigrid to be operated on.
void DeleteNumb(IdxDescCL& idx, MultiGridCL& MG);

} // end of namespace DROPS


#endif

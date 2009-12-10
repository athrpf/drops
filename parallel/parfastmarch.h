/// \file parfastmarch.h
/// \brief parallel version of fast marching method (line by line parallelization)
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier

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

/// This method is, from a performance point of view, not good. So, this
/// file may become deprecated.

#ifndef DROPS_PARFASTMARCH_H
#define DROPS_PARFASTMARCH_H

#include "parallel/parallel.h"
#include "parallel/exchange.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/spmat.h"
#include <vector>
#include <set>
#include <list>

namespace DROPS
{

typedef std::pair<bool,IdxT> FMIdxT;                    ///< Coupling of flag if DoF is a ghost and index

//------------------------
// R E P R  T E T R A  C L
//------------------------
/// \brief Store a representative tetra
class ReprTetraCL
/// A teraheader is stored by the four DoF. Each DoF is expanded by a flag, if the
/// DoF is a ghost or a local stored DoF
{
  private:
    SArrayCL<FMIdxT,4> data_;                           // Store a tetra by DoF

  public:
    ReprTetraCL() : data_(FMIdxT(false,NoIdx)) {}

    inline IdxT   operator[] (Uint) const;              ///< Get DoF
    inline bool   IsGhost    (Uint) const;              ///< Ask if DoF is a ghost
    inline void   Set        (Uint,const FMIdxT&);      ///< Set DoF and ghost-flag
    inline FMIdxT Get        (Uint) const;              ///< Get DoF and ghost-flag
           bool   HasGhost   () const;                  ///< Is any DoF a ghost
    /// Check if the DoF are the same (no distinction of ghost)
    friend bool operator==(const ReprTetraCL&, const ReprTetraCL&);
};

typedef std::pair<int,IdxT>  GhostIdxT;                 ///< Coupling of proc id and index on other proc
typedef std::list<GhostIdxT> GhostsList;                ///< List of all ghosts of a DoF
typedef std::vector<ReprTetraCL> VertexNeighT;          ///< All neighbor tetras of a single DoF
enum FMMark { Finished= 0, Close= 1, Far= 2};           ///< Types for verts (handled, next to handled, far away)


//--------------------------
// F M  T R A N S F E R  C L
//--------------------------
/// \brief Class for transfering within the fast marching algorithm
template<typename ExCL>
class FMTransferCL
/// This class initialize the parallel data structure before the fastmarching algortihm and handles the transfers
/// within the fastmarching algorithm.
{
  public:
    typedef std::list<ReprTetraCL> XferTetraList;           ///< List of tetras that shoul be transfered
    typedef std::list<IdxT>        IdxListT;                ///< List of hasghost DoF, that should be updated
    typedef std::set<int>          ProcSet;                 ///< Set of procs

    struct CoupFlagIdxS{                                    ///< Coupling of ghost-flag and index for transfer a ghost tetra or idx and type for updating type
        char flag;                                          ///< flag (ghost or type)
        IdxT dof;                                           ///< index
    };
    struct CoupCoordValS{                                   ///< Coupling of coordinate and levelset value for init ghosts
        double coord[3];                                    ///< coordinate
        double val;                                         ///< levelset value
        byte   typ;                                         ///< type of ghost
    };
    struct CoupIdxValS{                                     ///< Coupling of ghost index and levelset value for ghost updates
        char flag;                                          ///< flag if ghost or not
        IdxT dof;                                           ///< ghost index
        double val;                                         ///< levelset value
    };
    struct CoupIdxTypeS{                                    ///< Coupling of index and type
        char flag;                                          ///< flag if ghost or not
        IdxT dof;                                           ///< ghost index
        char val;                                           ///< type of DoF
    };
    const IdxT StopIdx;                                     ///< Stop index (NoIdx-1)

  private:
    ProcCL::DatatypeT           FlagIdxMPIT_;               // MPI Datatype for CoupFlagIdxS
    ProcCL::DatatypeT           CoordValMPIT_;              // MPI Datatype for CoupCoordValS (Can be also implemented by MPI::Datatye::Get_contiguous, but I find this looks better)
    ProcCL::DatatypeT           IdxValMPIT_;                // MPI Datatype for CoupIdxValS
    ProcCL::DatatypeT           IdxTypeMPIT_;               // MPI Datatype for CoupIdxTypeS
    VecDescCL&                  v_;                         // values for all local DoF
    VectorBaseCL<Point3DCL>&    Coord_;                     // coordinates of all local DoF
    VectorBaseCL<VertexNeighT>& neigh_;                     // tetras that owns a vert
    VectorBaseCL<byte>&         Typ_;                       // Each vert is Finished, Close or Far
    std::vector<double>&        DataExt_;                   // datas on ghost DoF
    std::vector<Point3DCL>&     CoordExt_;                  // coordinates of ghost DoF
    std::vector<byte>&          TypExt_;                    // types on ghost DoF
    VectorBaseCL<GhostsList>&   HasGhost_;                  // DoFs that has ghost on other procs
    ExCL&                       ex_;                        // ExchangeCL for computing external indices and distribution of DoF
    XferTetraList               ToXferTetra_;               // List of all tetras to be send
    IdxListT                    ToUpdateGhost_;             // List of hasghost DoF that has be changed
    IdxListT                    ToUpdateType_;              // List of updated types

    typename ExCL::ProcNumCT    neighs_;                    // neighbor procs
    int                         tag_;                       // Tag used inside this class (set to 3001)

    void CreateFlagIdxMPIT();                               // Create the MPI Datatypes
    void CreateCoordValMPIT();
    void CreateIdxValMPIT();
    void CreateIdxTypeMPIT();

    bool    HasDist(const ReprTetraCL&) const;              // Check if any DoF in RepTetra is distributed
    ProcSet GetProcs(const ReprTetraCL&) const;             // Get Procs that have a common DoF of a tetra

  public:
    FMTransferCL(VecDescCL&, VectorBaseCL<Point3DCL>&,
                 VectorBaseCL<VertexNeighT>&, VectorBaseCL<byte>&,
                 std::vector<double>&, std::vector<Point3DCL>&,
                 std::vector<byte>&,
                 VectorBaseCL<GhostsList>&, ExCL&);
    ~FMTransferCL();

    void MarkTetraForXfer(const ReprTetraCL&);              ///< Mark a tetra for transfer
    void XferTetras();                                      ///< Send and recieve all tetras that are marked

    void ChangedHasGhost(IdxT);                             ///< Remember changed HasGhost-DoF
    void UpdateGhosts();                                    ///< Update all changed values on HasGhost-DoFs

    void MarkChangeType(IdxT);                              ///< Mark DoF if type has changed
    void UpdateType();                                      ///< Update all type changings

    IdxT FindTrial(IdxT) const;                             ///< Find DoF with minimal value within Close
};


//------------------------------
// P A R  F A S T M A R C H  C L
//------------------------------
/// \brief Class for performing the parallel fastmarching algorithm
template<typename ExCL>
class ParFastMarchCL
{
  private:
    MultiGridCL&               MG_;                     // multigrid on which the fast marching method should be applied
    ExCL&                      ex_;                     // ExchangeCL for computing external indices and distribution of DoF
    const IdxT                 size_;                   // Number of unknowns of the levelset function

    // local values
    VecDescCL&                 v_;                      // Reference on the levelset VecDescCL
    std::set<IdxT>             Close_;                  // Verts that have finished neighbors
    VectorBaseCL<byte>         Typ_;                    // Each vert is Finished, Close or Far
    VectorBaseCL<Point3DCL>    Coord_;                  // coordinates of all verts
    VectorBaseCL<VertexNeighT> neigh_;                  // tetras that owns a vert
    VectorBaseCL<GhostsList>   HasGhost_;               // has another proc a ghost copy of this DoF, then store the representative tetras
    VectorCL                   Old_;                    // old values befor reparametrization

    // distributed values
    std::vector<Point3DCL>     CoordExt_;               // coordinates of ghost DoF
    std::vector<double>        DataExt_;                // datas on ghost DoF
    std::vector<byte>          TypExt_;                 // types on ghost DoF

    FMTransferCL<ExCL>         FMTransfer_;             // Transfering for fastmarching algorithm
    bool initpar_;                                      // is the parallel structure build

    void      UpdateHasGhost(IdxT);                     // Mark has ghost index, that the value is changed
    bool      HasGhost(IdxT) const;                     // Check if a DoF has a ghost DoF

    Point3DCL GetCoord(const FMIdxT&) const;            // Get coordinate of a DoF
    double    GetData(const FMIdxT&) const;             // Get levelset value of a DoF
    void      SetData(double, IdxT);                    // Set levelset value of a DoF
    byte      GetType(const FMIdxT&) const;             // Get type of a DoF
    void      SetType(FMMark, IdxT);                    // Set Mark on a dof

    double CompValueProj( IdxT, int, const FMIdxT upd[3]) const;
    void   Update(const FMIdxT&);                       // Update levelset value on a DoF
    IdxT   FindTrial() const;                           // Find next DoF
    void   UpdateGhost();                               // Update all Ghosts

    // Debugging Purpose
    bool CheckHasGhosts();                              // Check if a proc has a Ghost than the proc must not have this DoF as dist too
    bool AllPos();                                      // Check if all values are greater or equal zero

  public:
    ParFastMarchCL( MultiGridCL&, VecDescCL&, ExCL&);

    void InitCoord();                                   ///< Init coordinates
    void InitNeighAndPar();                             ///< Init neighborship for all DoF and parallel structure
    void InitZero( bool ModifyZero);                    ///< Init tetras around zero-level of the levelset function
    void InitClose();                                   ///< Init the "close-set"
    void Reparam( bool ModifyZero= true);               ///< Reparametrization of the levelset function
    void RestoreSigns();                                ///< Restore old signs
};

template<typename ConT>
  inline bool IsIn(const ConT &T, const typename ConT::value_type);

} // end of namespace DROPS

// File, where the inline an template-functions are declared
#include "parallel/parfastmarch.tpp"

#endif

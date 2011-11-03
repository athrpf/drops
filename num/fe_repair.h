/// \file fe_repair.h
/// \brief Repair-classes with respect to multigrid-changes for finite-elements.
/// \author LNM RWTH Aachen: Joerg Grande; SC RWTH Aachen:

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

/// \todo Remove RepairAfterRefineP2 from num/fe.*
/// \todo Do the same for the P1-case.

#ifndef DROPS_FE_REPAIR_H
#define DROPS_FE_REPAIR_H

#include <tr1/unordered_map>
#include <tr1/unordered_set>

#include "misc/container.h"
#include "geom/topo.h"
#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/discretize.h"

namespace DROPS
{

///\brief A vector of dof with their barycentric position with respect to the parent-tetra of the LocalP2-functions in RepairP2DataCL::data. The first pair-component is the position in the new VecDescCL.
typedef std::vector<std::pair<size_t,BaryCoordCL> > AugmentedDofVecT;

/// \brief Decide whether a barycentric point is in the closed reference tetra.
/// Optionally, a component-wise error can be specified.
inline bool
contained_in_reference_tetra (const BaryCoordCL& p, double eps= 0.)
{
    return p[0] >= -eps  && p[1] >= -eps && p[2] >= -eps && p[3] >= -eps;
}

/// \brief Repair-data for a P2-FE-function on a tetra, cf. RepairP2CL.
template <class ValueT>
struct RepairP2DataCL
{
    typedef ValueT value_type;
    typedef std::pair<Ubyte, LocalP2CL<value_type> > ChildDataT;
    typedef std::vector<ChildDataT> ChildVecT;

    ChildVecT data; ///< tuple of (old child number per topo.h, corresponding numerical P2-data)

    /// \brief Evaluate all augmented dof in dof on the repair-data in data and put the values into newdata.
    void repair (AugmentedDofVecT& dof, VectorCL& newdata) const;
};

/// \brief Repair a P2-FE-function after refinement.
///
/// The repair works as follows:
/// Before the refinement-algo, but after all marks for the refinement-algo have
/// been set, the data from all tetras, which could be deleted, are saved in
/// RepairP2Data-structs. To access these, the address of the parent is used,
/// because the refinement-algo does not remove a parent.
/// 
/// Tetras, which could be removed are
///     1) those with a mark for removement and level > 0,
///     2) irregular leaves (which always have level > 0).
/// 
/// There are 4 ways in which the refinement-algo can modify a tetra t in
/// triangulation l. The helpers of repair() below correspond correspond to these.
///     1) Deletion: t is removed and its parent p is now in l. The repair-data on p
/// (which is now in l) identifies this situation.
/// 
///     2) No change: t remains in l. This comes in two flavors: a) t is in lvl 0 and
/// there is no repair-data for l. b) t is in a higher level and there is no
/// repair-data on the parent and the grand-parent. (b) is determined as default
/// after ruling out (1), (2a), (3), (4).
/// 
///     3) Changed refinement of parent: t is deleted, but p is refined differently
/// again. The repair-data on the parent identifies this situation.
/// 
///     4) Refinement of t (or its regular replacement, if t was irregular), where the
/// children c are in l: a) If there is a grand-parent gp, this identifies the
/// situation. Its repair-data is used. b) t is in grid-level 1. There is no
/// grand-parent to decide, whether t was just created or unchanged as in case (2).
/// The tie is broken by recording the leave-tetras in level 0 in pre_refine(). If
/// the parent of t is such a leaf, then t is newly created, otherwise it remained
/// unchanged.
///
/// The algorithm is optimal in the sense, that the repaired data is always a
/// quadratic interpolant on the new triangulation of the original data (in contrast
/// to the earlier approaches used in Drops). Further, in cases (2) and (4), the
/// original function is interpolated exactly.
///
/// In principle, this class could use the accumulator-pattern (e.g., as base-class).
/// As this slightly complicates its use and because repair is performed seldom,
/// we leave this as future work.
template <class ValueT>
class RepairP2CL
{
  public:
    typedef ValueT value_type;

  private:
    typedef std::tr1::unordered_map<const TetraCL*, RepairP2DataCL<value_type> > RepairMapT;
    typedef std::tr1::unordered_set<const TetraCL*>                              TetraSetT;

    RepairMapT parent_data_;
    TetraSetT  level0_leaves_;

    BaryCoordCL p2_dof_[10]; ///< The bary-coordinates of the P2-dof.

    const MultiGridCL& mg_;            ///< Multigrid to operate on
    const VecDescCL& old_vd_;          ///< original data to be repaired
    const BndDataCL<value_type>& bnd_; ///< boundary-data for the old and new VecDescCL

    VecDescCL* new_vd_;                ///< VecDescCL to be repaired; set in repair().
    std::vector<bool> repair_needed_;   ///< Memoize dof that must be repaired. Used only in repair() and its helpers.

    bool repair_needed    (size_t dof) { return repair_needed_[dof]; }
    void mark_as_repaired (size_t dof) { repair_needed_[dof]= false; }

    AugmentedDofVecT collect_unrepaired_dofs (const TetraCL& t); ///< collect dofs with repair_needed().
    void unchanged_refinement    (const TetraCL& t); ///< use data from t for copying
    void regular_leaf_refinement (const TetraCL& t); ///< use data from t for repair
    void unrefinement            (const TetraCL& t, const RepairP2DataCL<ValueT>& t_data);  ///< use repair-data from the tetra itself
    void changed_refinement      (const TetraCL& t, const RepairP2DataCL<ValueT>& p_data);  ///< use repair-data from the parent
    void genuine_refinement      (const TetraCL& t, const RepairP2DataCL<ValueT>& gp_data); ///< use repair-data from the grand-parent

  public:
    /// \brief Initializes the data to be repaired on mg and calls pre_refine().
    RepairP2CL (const MultiGridCL& mg, const VecDescCL& old, const BndDataCL<value_type>& bnd);

    /// \brief Saves data from possibly deleted tetras in parent_data_ and level0_leaves_.
    void pre_refine ();

    /// \brief Repair old_vd with the help of the saved data and store the result in new_vd.
    /// new_vd is assumed to contain a valid numbering and a vector of corresponding size.
    /// The new numbering must use the same boundary-data as the old one.
    void repair (VecDescCL& new_vd);
};

} // end of namespace DROPS

#include "num/fe_repair.tpp"

#endif

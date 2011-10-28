/// \file fe_repair.tpp
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

namespace DROPS
{

/// RepairP2DataCL

template <class ValueT>
void RepairP2DataCL<ValueT>::repair (AugmentedDofVecT& dof, VectorCL& newdata) const
{
    AugmentedDofVecT::iterator d;
    BaryCoordCL tmp( Uninitialized);
    for (typename ChildVecT::const_iterator old_ch= data.begin(); !(dof.empty() || old_ch == data.end()); ++old_ch) {
        const SMatrixCL<4,4>& to_old_child= parent_to_child_bary( old_ch->first);
        d= dof.begin();
        while (d != dof.end()) {
            tmp= to_old_child*d->second;
            if (contained_in_reftetra( tmp)) {
                DoFHelperCL<value_type, VectorCL>::set( newdata, d->first, old_ch->second( tmp));
                d= dof.erase( d);
            }
            else
                ++d;
        }
    }
    if (!dof.empty())
        throw DROPSErrCL("RepairP2DataCL::repair: Could not locate all new dof.\n");
}


/// RepairP2CL

template <class ValueT>
const double RepairP2CL<ValueT>::UnrepairedDofC= std::numeric_limits<double>::quiet_NaN();

template <class ValueT>
  RepairP2CL<ValueT>::RepairP2CL (const MultiGridCL& mg, const VecDescCL& old, const BndDataCL<value_type>& bnd)
        : mg_( mg), old_vd_ ( old), bnd_( bnd)
{
    for (Uint i= 0; i < NumVertsC; ++i)
        p2_dof_[i]= std_basis<4>( i + 1);
    for (Uint i= 0; i < NumEdgesC; ++i)
        p2_dof_[i + NumVertsC]= 0.5*(std_basis<4>( VertOfEdge( i, 0) + 1) + std_basis<4>( VertOfEdge( i, 1) + 1));

    pre_refine();
}

template <class ValueT>
  void
  RepairP2CL<ValueT>::pre_refine ()
{
    Uint lvl= old_vd_.RowIdx->TriangLevel();
    LocalP2CL<value_type> lp2;
    DROPS_FOR_TRIANG_CONST_TETRA( mg_, lvl, it) {
        if (!it->IsUnrefined())
            continue;
        if (it->GetLevel() > 0) {
            // These could be deleted by the refinement algo
            /// \todo To store less, one can generally add "if ( it->IsMarkedForRemovement() || !it->IsRegular())" here. However, the TetrabuilderCL::BogoReMark-algorithm does not obey the user-marks, thus we postpone this optimization.
            const TetraCL* p= it->GetParent();
            const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &*it) - p->GetChildBegin();
            const RefRuleCL& rule= p->GetRefData();
            lp2.assign_on_tetra( *it, old_vd_, bnd_);
            parent_data_[p].data.push_back( std::make_pair( rule.Children[ch], lp2));
        }
        else {
            // Leaves in level 0 can give birth to children, which cannot easily be distinguished from tetras, which just remained over one refinement step, cf. case (2) and (4b) in repair(). We memoize the leaves in level 0, as they tend to be few.
            level0_leaves_.insert( &*it);
        }
    }
}

template <class ValueT>
  AugmentedDofVecT
  RepairP2CL<ValueT>::collect_unrepaired_dofs (const TetraCL& t)
{
    LocalNumbP2CL n_new( t, *new_vd_->RowIdx);
    AugmentedDofVecT dof;
    dof.reserve( 10);
    for (Uint i= 0; i < 10; ++i)
        if (n_new.WithUnknowns( i) && repair_needed( new_vd_->Data[n_new.num[i]]))
            dof.push_back( std::make_pair( n_new.num[i], p2_dof_[i]));
    return dof;
}

template <class ValueT>
void RepairP2CL<ValueT>::unchanged_refinement (const TetraCL& t)
{
    const VectorCL& olddata= old_vd_.Data;
          VectorCL& newdata= new_vd_->Data;
    LocalNumbP2CL n_old( t, *old_vd_.RowIdx);
    LocalNumbP2CL n_new( t, *new_vd_->RowIdx);
    for (Uint i= 0; i < 10; ++i)
        if (n_new.WithUnknowns( i) && repair_needed( newdata[n_new.num[i]])) {
            Assert( n_old.WithUnknowns( i), DROPSErrCL( "TetraRepairP2CL::unchanged_refinement: "
                "Old and new function must use the same boundary-data-types.\n"), DebugNumericC);
            const value_type& tmp= DoFHelperCL<value_type, VectorCL>::get( olddata, n_old.num[i]);
            DoFHelperCL<value_type, VectorCL>::set( newdata, n_new.num[i], tmp);
        }
}

template <class ValueT>
void RepairP2CL<ValueT>::regular_leaf_refinement (const TetraCL& t)
{
    const AugmentedDofVecT& dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    const SMatrixCL<4,4>& T= child_to_parent_bary( p->GetRefData().Children[ch]);

    LocalP2CL<value_type> oldp2;
    oldp2.assign_on_tetra( *p, old_vd_, bnd_);
    for (AugmentedDofVecT::const_iterator d= dof.begin(); d != dof.end(); ++d)
        DoFHelperCL<value_type, VectorCL>::set( new_vd_->Data, d->first, oldp2( T*d->second));
}

template <class ValueT>
void RepairP2CL<ValueT>::genuine_refinement (const TetraCL& t, const RepairP2DataCL<ValueT>& repairdata)
{
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    const SMatrixCL<4,4>& T= child_to_parent_bary( p->GetRefData().Children[ch]);
    const TetraCL* gp= p->GetParent();
    const Ubyte gpch= std::find( gp->GetChildBegin(), gp->GetChildEnd(), p) - gp->GetChildBegin();
    const SMatrixCL<4,4>& S= child_to_parent_bary( gp->GetRefData().Children[gpch]);
    for (AugmentedDofVecT::iterator d= dof.begin(); d != dof.end(); ++d)
        d->second= S*(T*d->second);
    repairdata.repair( dof, new_vd_->Data);
}

template <class ValueT>
void RepairP2CL<ValueT>::unrefinement (const TetraCL& t, const RepairP2DataCL<ValueT>& repairdata)
{
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    QRDecompCL<4> T;
    BaryCoordCL tmp;
    typedef typename RepairP2DataCL<value_type>::ChildVecT ChildVecT;
    repairdata.repair( dof, new_vd_->Data);
}

template <class ValueT>
void RepairP2CL<ValueT>::changed_refinement (const TetraCL& t, const RepairP2DataCL<ValueT>& repairdata)
{
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    const SMatrixCL<4,4>& to_parent= child_to_parent_bary( p->GetRefData().Children[ch]);
    for (AugmentedDofVecT::iterator d= dof.begin(); d != dof.end(); ++d)
        d->second= to_parent*d->second;
    repairdata.repair( dof, new_vd_->Data);
}

template <class ValueT>
  void
  RepairP2CL<ValueT>::repair (VecDescCL& new_vd)
{
    new_vd_= &new_vd;

    const Uint lvl= new_vd_->RowIdx->TriangLevel();
    Assert( lvl == old_vd_.RowIdx->TriangLevel() || lvl ==  old_vd_.RowIdx->TriangLevel() - 1,
        DROPSErrCL( "RepairP2CL<ValueT>::repair: Different levels\n"), DebugNumericC);
    if (lvl == old_vd_.RowIdx->TriangLevel() - 1)
        std::cout << "old level: " << old_vd_.RowIdx->TriangLevel() << " mg_.GetLastLevel(): " << mg_.GetLastLevel() << '\n';

    VectorCL& newdata= new_vd_->Data;
    newdata= UnrepairedDofC;

    DROPS_FOR_TRIANG_CONST_TETRA( mg_, lvl, t) {
        if (parent_data_.count( &*t) == 1) // Case 1
            unrefinement( *t, parent_data_[&*t]);
        // From here on, t has no parent-data itself.
        else if (t->GetLevel() == 0) // Case 2
            // If t had arrived in the triangulation via coarsening, it would have had parent-data.
            unchanged_refinement( *t);
        else { // t has no parent-data and t->GetLevel() > 0
            const TetraCL* p= t->GetParent();
            if (parent_data_.count( p) == 1) // Case 3
                changed_refinement( *t, parent_data_[p]);
            else { // t has no repair-data, t->GetLevel() > 0, and p has no repair-data
                if ((p->GetLevel() > 0 && parent_data_.count( p->GetParent()) == 1))
                    genuine_refinement( *t, parent_data_[p->GetParent()]); // Case (4a).
                else if (level0_leaves_.count( p) == 1)
                    regular_leaf_refinement( *t); // Case (4b).
                else // Case (2)
                    unchanged_refinement( *t);
            }
        }
    }
}

} // end of namespace DROPS

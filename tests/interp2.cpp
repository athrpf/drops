/// \file interp2.cpp
/// \brief tests implementation of RepareAfterRefineP2
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "num/discretize.h"
#include <tr1/unordered_set>


namespace DROPS {

bool
contained_in_reftetra( const BaryCoordCL& p, double eps= 0.)
{
    return p[0] >= -eps  && p[1] >= -eps && p[2] >= -eps && p[3] >= -eps;
}

void
to_parent_bary (Uint ch, SMatrixCL<4,4>& T)
{
    const ChildDataCL& child= GetChildData( ch);
    for (Uint i= 0; i < NumVertsC; ++i) {
        const Uint v= child.Vertices[i];
        if (v < NumVertsC)
            T.col( i, std_basis<4>( v + 1));
        else
            T.col( i, 0.5*(std_basis<4>( VertOfEdge( v - 4, 0) + 1) + std_basis<4>( VertOfEdge( v - 4, 1) + 1)));
    }
}

void
to_child_bary (Uint ch, QRDecompCL<4>& T)
{
    to_parent_bary( ch, T.GetMatrix());
    T.prepare_solve();
}

template <class ValueT>
struct RepairP2DataCL
{
    typedef ValueT value_type;
    typedef std::pair<Ubyte, LocalP2CL<value_type> > ChildDataT;
    typedef std::vector<ChildDataT> ChildVecT;

    ChildVecT data; ///< tuple of (old child number per topo.h, corresponding numerical P2-data)

    typedef std::vector<std::pair<size_t,BaryCoordCL> > AugmentedDofVecT;
    void repair (AugmentedDofVecT& dof, VectorCL& newdata) const;
};

template <class ValueT>
void RepairP2DataCL<ValueT>::repair (AugmentedDofVecT& dof, VectorCL& newdata) const
{
    QRDecompCL<4> to_old_child;
    BaryCoordCL tmp;
    for (typename ChildVecT::const_iterator old_ch= data.begin(); !(dof.empty() || old_ch == data.end()); ++old_ch) {
        to_child_bary( old_ch->first, to_old_child);
        AugmentedDofVecT::iterator d= dof.begin();
        while (d != dof.end()) {
            tmp= d->second;
            to_old_child.Solve( tmp);
            if (contained_in_reftetra( tmp, 4.*std::numeric_limits<double>::epsilon())) {
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

bool is_repaired (double d) { return !std::isnan( d); }

const double UnrepairedDofC= std::numeric_limits<double>::quiet_NaN();

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
/// after ausschliessen (1), (2a), (3), (4).
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
template <class ValueT>
class RepairP2CL
{
  public:
    typedef ValueT value_type;

  private:
    typedef std::tr1::unordered_map<const TetraCL*, RepairP2DataCL<value_type> > RepairMapT;
    typedef std::tr1::unordered_set<const TetraCL*>                              TetraSetT;
    typedef std::vector<std::pair<size_t,BaryCoordCL> >                          AugmentedDofVecT;

    RepairMapT parent_data_;
    TetraSetT  level0_leaves_;

    BaryCoordCL p2_dof_[10];

    const MultiGridCL& mg_;
    const VecDescCL& old_vd_;
    const BndDataCL<value_type>& bnd_;

    VecDescCL* new_vd_;

    void pre_refine ();

    AugmentedDofVecT collect_unrepaired_dofs (const TetraCL& t); ///< collect dofs with !is_repaired.
    void unchanged_refinement    (const TetraCL& t); ///< use data from t for copying
    void regular_leaf_refinement (const TetraCL& t); ///< use data from t for repair
    void unrefinement            (const TetraCL& t, const RepairP2DataCL<ValueT>& t_data);  ///< use repair-data from the tetra itself
    void changed_refinement      (const TetraCL& t, const RepairP2DataCL<ValueT>& p_data);  ///< use repair-data from the parent
    void genuine_refinement      (const TetraCL& t, const RepairP2DataCL<ValueT>& gp_data); ///< use repair-data from the grand-parent

  public:
    RepairP2CL (const MultiGridCL& mg, const VecDescCL& old, const BndDataCL<value_type>& bnd);

    void repair (VecDescCL& new_vd);
};

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
            /// \todo To store less, one can add " && (it->IsMarkedForRemovement() || !it->IsRegular())"
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

template <class ValueT>
  typename RepairP2CL<ValueT>::AugmentedDofVecT
  RepairP2CL<ValueT>::collect_unrepaired_dofs (const TetraCL& t)
{
    LocalNumbP2CL n_new( t, *new_vd_->RowIdx);
    AugmentedDofVecT dof;
    dof.reserve( 10);
    for (Uint i= 0; i < 10; ++i)
        if (n_new.WithUnknowns( i) && !is_repaired( new_vd_->Data[n_new.num[i]]))
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
        if (n_new.WithUnknowns( i) && !is_repaired( newdata[n_new.num[i]])) {
            Assert( n_old.WithUnknowns( i), DROPSErrCL( "TetraRepairP2CL::unchanged_refinement: "
                "Old and new function must use the same boundary-data-types.\n"), DebugNumericC);
            const value_type& tmp= DoFHelperCL<value_type, VectorCL>::get( olddata, n_old.num[i]);
            DoFHelperCL<value_type, VectorCL>::set( newdata, n_new.num[i], tmp);
        }
}

template <class ValueT>
void RepairP2CL<ValueT>::regular_leaf_refinement (const TetraCL& t)
{
    typedef std::vector<std::pair<size_t,BaryCoordCL> > AugmentedDofVecT;
    const AugmentedDofVecT& dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    SMatrixCL<4,4> T( Uninitialized);
    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    to_parent_bary( p->GetRefData().Children[ch], T);

    LocalP2CL<value_type> oldp2;
    oldp2.assign_on_tetra( *p, old_vd_, bnd_);
    for (AugmentedDofVecT::const_iterator d= dof.begin(); d != dof.end(); ++d)
        DoFHelperCL<value_type, VectorCL>::set( new_vd_->Data, d->first, oldp2( T*d->second));
}

template <class ValueT>
void RepairP2CL<ValueT>::genuine_refinement (const TetraCL& t, const RepairP2DataCL<ValueT>& repairdata)
{
    typedef std::vector<std::pair<size_t,BaryCoordCL> > AugmentedDofVecT;
    AugmentedDofVecT dof= collect_unrepaired_dofs( t);
    if (dof.empty())
        return;

    SMatrixCL<4,4> T( Uninitialized);
    const TetraCL* p= t.GetParent();
    const Ubyte ch= std::find( p->GetChildBegin(), p->GetChildEnd(), &t) - p->GetChildBegin();
    to_parent_bary( p->GetRefData().Children[ch], T);
    const TetraCL* gp= p->GetParent();
    const Ubyte gpch= std::find( gp->GetChildBegin(), gp->GetChildEnd(), p) - gp->GetChildBegin();
    SMatrixCL<4,4> S( Uninitialized);
    to_parent_bary( gp->GetRefData().Children[gpch], S);
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
    SMatrixCL<4,4> to_parent( Uninitialized);
    to_parent_bary( p->GetRefData().Children[ch], to_parent);
    for (AugmentedDofVecT::iterator d= dof.begin(); d != dof.end(); ++d)
        d->second= to_parent*d->second;
    repairdata.repair( dof, new_vd_->Data);
}

} // end of namespace DROPS



using namespace DROPS;

typedef double (*fun_ptr)(const SVectorCL<3>&);

enum  OutputModeT { SILENT, NOISY };


double f(const SVectorCL<3>& p)
{ return p[0]*p[0] +10.*p[1]*p[1] +100.*p[2]*p[2] +1000.*p[0]*p[1] +10000.*p[0]*p[2] +100000.*p[1]*p[2]; }

double g(const SVectorCL<3>& p)
{  return p[0] +10.*p[1] +100.*p[2]+1000.; }

double h(const SVectorCL<3>& p)
{  return std::sin(M_PI*p[0])*std::sin(M_PI*p[1])*std::sin(M_PI*p[2]); }

double g2(const DROPS::Point3DCL& p)
{
    return (-1.)*p[0];
}

void MarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop(DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte( 0.5);

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It( mg.GetTriangTetraBegin(maxLevel)),
             ItEnd( mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It) {
        if ( (GetBaryCenter( *It)-Mitte).norm() <= std::max( 0.2, 1.5*std::pow( It->GetVolume(), 1.0/3.0)) )
            It->SetRemoveMark();
    }
}

typedef NoBndDataCL<double> BndCL;
BndCL Bnd;

void SetFun(VecDescBaseCL<VectorCL>& vd, MultiGridCL& mg, fun_ptr f)
{
    vd.Data.resize( vd.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun( &vd, &Bnd, &mg);
    const Uint lvl= vd.GetLevel();
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(lvl),
         theend= mg.GetTriangVertexEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( sit->GetCoord()));
    }
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(lvl),
         theend= mg.GetTriangEdgeEnd(lvl); sit!=theend; ++sit) {
        fun.SetDoF( *sit, f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5));
    }
}


int CheckResult(DROPS::P2EvalCL<double, BndCL,
                const DROPS::VecDescCL>& fun, fun_ptr f, OutputModeT om, double eps= 0.)
{
    int ret= 0;
    const VertexCL* v= 0;
    const EdgeCL* e= 0;
    const DROPS::MultiGridCL& mg= fun.GetMG();
    const DROPS::Uint trilevel= fun.GetLevel();
    if (om!=SILENT) std::cout << "Verts:" << std::endl;
    double diff, emaxdiff= 0., vmaxdiff= 0., vval= 1., eval;
    for (MultiGridCL::const_TriangVertexIteratorCL sit=mg.GetTriangVertexBegin( trilevel),
         theend= mg.GetTriangVertexEnd( trilevel); sit!=theend; ++sit) {
        diff= fun.val( *sit) - f( sit->GetCoord());
        if (std::abs( diff) > eps && std::abs(diff) > vmaxdiff) {
            ++ret;
            vmaxdiff= std::abs(diff);
            v= &*sit;
            vval= fun.val( *sit);
        }
    }
    if (om!=SILENT) std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::const_TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin( trilevel),
         theend= mg.GetTriangEdgeEnd( trilevel); sit!=theend; ++sit) {
        diff = fun.val( *sit, .5) - f( (sit->GetVertex( 0)->GetCoord() + sit->GetVertex( 1)->GetCoord())*0.5);
        if (std::abs( diff) > eps && std::abs(diff) > emaxdiff) {
            ++ret;
            emaxdiff= std::abs(diff);
            e= &*sit;
            eval= fun.val( *sit);
        }
    }
    if (om!=SILENT) {
        if (vval == 0.) vval= 1.;
        std::cout << "maximale Differenz Vertices: " << vmaxdiff << " (exakter Wert: " << vval << ") auf\n";
        if (v) v->DebugInfo( std::cout);
        if (eval == 0.) eval= 1.;
        std::cout << "maximale Differenz Edges: " << emaxdiff << " (exakter Wert: " << eval << ") auf\n";
        if (e) e->DebugInfo( std::cout);
        std::cout << std::endl;
    }
    return ret;
}


DROPS::Uint Rule(DROPS::Uint r)
{
    return r < 64 ? r : 127;
}

// True, iff every bit in p is also set in q.
bool SubSuperPattern(DROPS::Uint p, DROPS::Uint q)
{
    return !((~q) & p);
}


// Checks every possible tetra-modification.
int TestReMark()
{
    BndDataCL<> bnd( 4);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair on single tetra-combinations:\n";
    int ttt, ret= 0;
    for (DROPS::Uint i= 0; i<=64; ++i) {
        for (DROPS::Uint j= 0; j<=64; ++j) {
            DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
            // std::cout << Rule( i) << "\t-->\t" << Rule( j) << " ";
            DROPS::TetraBuilderCL tet( Rule( i), DROPS::std_basis<3>( 1),
                                                 DROPS::std_basis<3>( 2),
                                                 DROPS::std_basis<3>( 3),
                                                 DROPS::Point3DCL( 1.0));
            DROPS::MultiGridCL mg( tet);
            DROPS::IdxDescCL i0(P2_FE, Bnd), i1(P2_FE, Bnd);
            i0.CreateNumbering( mg.GetLastLevel(), mg);
            DROPS::VecDescCL v0, v1;
            v0.SetIdx( &i0);
            SetFun( v0, mg, f);
            // SetFun( v0, mg, g2);
            RepairP2CL<double> repairp2( mg, v0, bnd);
            tet.BogoReMark( mg, Rule( j));

            const Uint i1_Level= i0.TriangLevel() <= mg.GetLastLevel() ? i0.TriangLevel()
                                                                       : mg.GetLastLevel();
            i1.CreateNumbering( i1_Level, mg);
            v1.SetIdx( &i1);
            repairp2.repair( v1);
            DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
//            ttt= CheckResult( fun1, f, SILENT);
            ttt= CheckResult( fun1, f, SILENT, 1e-10);
            ret+= ttt;
            if (ttt != 0 && SubSuperPattern( Rule( i) & 63, Rule( j) & 63))
                std::cout << "Aerger: " << Rule( i) << "\t-->\t" << Rule( j) << " " << std::endl;
        }
    }
    return ret;
}


int TestRepairUniform()
{
    BndDataCL<> bnd( 6);
    int ret= 0;
    DROPS::BrickBuilderCL brick( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                                 DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                                 1, 1, 1);
    DROPS::MultiGridCL mg(brick);
    DROPS::IdxDescCL i0( P2_FE, Bnd), i1( P2_FE, Bnd);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for uniform refinement with quadratic function:\n";
    for (DROPS::Uint i=0; i<5; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx( &i0);
        SetFun( v0, mg, f);
        RepairP2CL<double> repairp2( mg, v0, bnd);
        MarkAll( mg);
        mg.Refine();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        v1.SetIdx( &i1);
        repairp2.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
    }

    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for uniform coarsening with quadratic function:\n";
    for (DROPS::Uint i=0; i<5; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg, f);
        RepairP2CL<double> repairp2( mg, v0, bnd);
        UnMarkAll( mg);
        mg.Refine();
        Uint i1_Level= i0.TriangLevel();
        if (mg.GetLastLevel() < i1_Level) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            --i1_Level;
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
        }
        i1.CreateNumbering( i1_Level, mg);
        v1.SetIdx( &i1);
        repairp2.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        if (mg.GetLastLevel() < i0.TriangLevel()) {
            Uint level= mg.GetLastLevel();
            DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllVertexBegin( level),
                                        mg.GetAllVertexEnd( level));
            DROPS::DeleteNumbOnSimplex( i0.GetIdx(), mg.GetAllEdgeBegin( level),
                                        mg.GetAllEdgeEnd( level));
        }
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllVertexBegin( i1.TriangLevel()),
                                    mg.GetAllVertexEnd( i1.TriangLevel()));
        DROPS::DeleteNumbOnSimplex( i1.GetIdx(), mg.GetAllEdgeBegin( i1.TriangLevel()),
                                    mg.GetAllEdgeEnd( i1.TriangLevel()));
    }
    return ret;
}

int TestRepair()
{
    BndDataCL<> bnd( 6);
    int ret= 0;
    DROPS::BrickBuilderCL brick( DROPS::std_basis<3>( 0), DROPS::std_basis<3>( 1),
                                 DROPS::std_basis<3>( 2), DROPS::std_basis<3>( 3),
                                 1, 1, 1);
    DROPS::IdCL<DROPS::VertexCL>::ResetCounter();
    DROPS::MultiGridCL mg(brick);
    DROPS::IdxDescCL i0( P2_FE, Bnd), i1( P2_FE, Bnd);
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for drop refinement with quadratic function:\n";
    for (DROPS::Uint i=0; i<8; ++i) {
        std::cout << "i: " << i;
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        std::cout << " i0.TriangLevel(): " << i0.TriangLevel();
        DROPS::VecDescCL v0, v1;
        v0.SetIdx( &i0);
        SetFun( v0, mg, f);
        RepairP2CL<double> repairp2( mg, v0, bnd);
        MarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        std::cout << " mg.GetLastLevel() after refine: " << mg.GetLastLevel();
        i1.CreateNumbering( i0.TriangLevel(), mg);
        std::cout << " i1.TriangLevel(): " << i1.TriangLevel();
        v1.SetIdx( &i1);
        repairp2.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
        std::cout << '\n';
    }

    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting repair for drop coarsening with quadratic function:\n";
    for (DROPS::Uint i=0; i<8; ++i) {
        i0.CreateNumbering( mg.GetLastLevel(), mg);
        DROPS::VecDescCL v0, v1;
        v0.SetIdx(&i0);
        SetFun(v0, mg, f);
        RepairP2CL<double> repairp2( mg, v0, bnd);
        UnMarkDrop( mg, mg.GetLastLevel());
        mg.Refine();
        Uint i1_Level= i0.TriangLevel();
        if (mg.GetLastLevel() < i1_Level) {
            std::cout << "Letztes Level entfernt!" << std::endl;
            --i1_Level;
        }
        else {
            std::cout << "Letztes Level behalten!" << std::endl;
        }
        i1.CreateNumbering( i1_Level, mg);
        v1.SetIdx( &i1);
        repairp2.repair( v1);
        DROPS::P2EvalCL<double, BndCL, const VecDescCL > fun1( &v1, &Bnd, &mg);
        ret+= CheckResult( fun1, f, NOISY, 1e-10);
        if (mg.GetLastLevel() < i0.TriangLevel())
            i0.DeleteNumbering( mg);
        i1.DeleteNumbering( mg);
    }
    return ret;
}


int TestInterpolateOld()
{
    std::cout << "\n-----------------------------------------------------------------"
                 "\nTesting interpolation:\n";
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 4, 4, 4);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 3, 3, 3);
    DROPS::MultiGridCL mg(brick);
    MarkDrop(mg, 0);
//    MarkAll( mg);
    mg.Refine();
    MarkDrop(mg, 1);
    mg.Refine();

    DROPS::IdxDescCL i0( DROPS::P2_FE, Bnd), i1( DROPS::P2_FE, Bnd);
    i0.CreateNumbering( 0, mg);
    i1.CreateNumbering( 1, mg);

    VecDescBaseCL<VectorCL> v0, v1;
    v0.SetIdx( &i0);
    v1.SetIdx( &i1);
    SetFun( v0, mg, f);

    v1.Data.resize( v1.RowIdx->NumUnknowns());
    P2EvalCL<double, BndCL, const VecDescBaseCL<VectorCL> > fun0( &v0, &Bnd, &mg);
    P2EvalCL<double, BndCL,VecDescBaseCL<VectorCL> > fun1( &v1, &Bnd, &mg);
    Interpolate( fun1, fun0);
//    DROPS::RepairAfterRefineP2( fun0, v1);
    std::cout << "Verts:" << std::endl;
    double diff;
    for (MultiGridCL::TriangVertexIteratorCL sit=mg.GetTriangVertexBegin(1),
         theend= mg.GetTriangVertexEnd(1); sit!=theend; ++sit) {
        diff= fun1.val(*sit) - f(sit->GetCoord());
        std::cout << diff << "\t";
        if (diff!=0.) return 1;
    }
    std::cout << "\n\nEdges:" << std::endl;
    for (MultiGridCL::TriangEdgeIteratorCL sit=mg.GetTriangEdgeBegin(1),
         theend= mg.GetTriangEdgeEnd(1); sit!=theend; ++sit) {
        diff = fun1.val( *sit, .5) - f( (sit->GetVertex(0)->GetCoord()+sit->GetVertex(1)->GetCoord())*.5);
        std::cout << diff << "\t";
        if (diff!=0.) return 1;
    }
    std::cout << std::endl;
    std::cout << std::endl << DROPS::SanityMGOutCL(mg) << std::endl;
    return 0;
}


int main ()
{
  try {
// Show the trafos and their inverses
// SMatrixCL<4,4> T( Uninitialized);
// QRDecompCL<4> S;
// SMatrixCL<4,4> Smat( Uninitialized);
// for (Uint i=0; i < NumAllChildrenC; ++i) {
//     to_parent_bary( i, T);
//     std::cout << i << '\n' << T << '\n';
//     to_child_bary( i, S);
//     for (Uint j= 0; j < 4; ++j) {
//         BaryCoordCL tmp= std_basis<4>( j + 1);
//         S.Solve( tmp);
//         Smat.col( j, tmp);
//     }
//     std::cout << i << '\n' << Smat << '\n';
// }
// return 0;

    int ret= TestRepairUniform();
    ret+= TestRepair();
    // ret+= TestInterpolateOld();
    return ret + TestReMark();
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

//**************************************************************************
// File:    adaptriang.tpp                                                 *
// Content: adaptive triangulation based on position of the interface      *
//          provided by the levelset function                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

namespace DROPS
{

template <class DistFctT>
void AdapTriangCL::MakeInitialTriang( DistFctT& Dist)
{
    TimerCL time;

    time.Reset();
    time.Start();
    const Uint min_ref_num= f_level_ - c_level_;
    Uint i;
    for (i=0; i<2*min_ref_num; ++i)
        ModifyGridStep( Dist);
    time.Stop();
    std::cout << "MakeTriang: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg_.GetLastLevel() << '\n';
    mg_.SizeInfo( std::cout);
}
    
template <class DistFctT>
bool AdapTriangCL::ModifyGridStep( DistFctT& Dist)
// One step of grid change; returns true if modifications were necessary,
// false, if nothing changed.
{
    bool modified= false;
    for (MultiGridCL::TriangTetraIteratorCL it= mg_.GetTriangTetraBegin(),
         end= mg_.GetTriangTetraEnd(); it!=end; ++it) 
    {
        double d= 1.;
        for (Uint j=0; j<4; ++j) 
            d= std::min( d, std::abs( GetValue( Dist, *it->GetVertex( j)) ));
        d= std::min( d, std::abs( GetValue( Dist, *it)));
        const Uint l= it->GetLevel();
	// In the shell:      level should be f_level_.
        // Outside the shell: level should be c_level_.
        const Uint soll_level= d<=width_ ? f_level_ : c_level_;

        if (l !=  soll_level)
	{ // tetra will be marked for refinement/remove
	    modified= true;
            if (l < soll_level) 
                it->SetRegRefMark();
            else // l > soll_level 
                it->SetRemoveMark();
        }
    }
    if (modified) 
	mg_.Refine();
    return modified;
}

template <class StokesT>
void AdapTriangCL::UpdateTriang( StokesT& NS, LevelsetP2CL& lset)
{
    TimerCL time;

    time.Reset();
    time.Start();
    VelVecDescCL  loc_v;
    VecDescCL     loc_p;
    VecDescCL     loc_l;
    VelVecDescCL *v1= &NS.v, 
                 *v2= &loc_v;
    VecDescCL    *p1= &NS.p,
                 *p2= &loc_p,
		 *l1= &lset.Phi,
		 *l2= &loc_l;
    IdxDescCL  loc_vidx( 3, 3), loc_pidx( 1), loc_lidx( 1, 1);
    IdxDescCL  *vidx1= v1->RowIdx,
               *vidx2= &loc_vidx,
	       *pidx1= p1->RowIdx,
	       *pidx2= &loc_pidx,
	       *lidx1= l1->RowIdx,
	       *lidx2= &loc_lidx;
    modified_= false;
    const Uint min_ref_num= f_level_ - c_level_;
    Uint i, LastLevel= mg_.GetLastLevel();

    for (i=0; i<2*min_ref_num; ++i)
    {            
	LevelsetP2CL::DiscSolCL sol= lset.GetSolution( *l1);
        if (!ModifyGridStep(sol))
            break;
        LastLevel= mg_.GetLastLevel();
        modified_= true;

        // Repair velocity
        std::swap( v2, v1);
        std::swap( vidx2, vidx1);
        NS.CreateNumberingVel( LastLevel, vidx1);
        if ( LastLevel != vidx2->TriangLevel) 
        {
            std::cout << "LastLevel: " << LastLevel
                      << " vidx2->TriangLevel: " << vidx2->TriangLevel << std::endl;
            throw DROPSErrCL( "AdapTriangCL::UpdateTriang: Sorry, not yet implemented.");
        }
        v1->SetIdx( vidx1);
        typename StokesT::DiscVelSolCL funvel= NS.GetVelSolution( *v2);
        RepairAfterRefineP2( funvel, *v1);
        v2->Clear();
        NS.DeleteNumberingVel( vidx2);

        // Repair pressure
        std::swap( p2, p1);
        std::swap( pidx2, pidx1);
        NS.CreateNumberingPr( LastLevel, pidx1);
        p1->SetIdx( pidx1);
        typename StokesT::DiscPrSolCL funpr= NS.GetPrSolution( *p2);
        RepairAfterRefineP1( funpr, *p1);
        p2->Clear();
        NS.DeleteNumberingPr( pidx2);

        // Repair levelset
        std::swap( l2, l1);
        std::swap( lidx2, lidx1);
        lset.CreateNumbering( LastLevel, lidx1);
        l1->SetIdx( lidx1);
        LevelsetP2CL::DiscSolCL funlset= lset.GetSolution( *l2);
        RepairAfterRefineP2( funlset, *l1);
        l2->Clear();
        lset.DeleteNumbering( lidx2);
    }
    // We want the solution to be in NS.v, NS.pr, lset.Phi
    if (v1 == &loc_v) 
    {
        NS.vel_idx.swap( loc_vidx);
        NS.pr_idx.swap( loc_pidx);
	lset.idx.swap( loc_lidx);
        NS.v.SetIdx( &NS.vel_idx);
        NS.p.SetIdx( &NS.pr_idx);
	lset.Phi.SetIdx( &lset.idx);

        NS.v.Data= loc_v.Data;
        NS.p.Data= loc_p.Data;
	lset.Phi.Data= loc_l.Data;
    }
    time.Stop();
    std::cout << "UpdateTriang: " << i
              << " refinements/interpolations in " << time.GetTime() << " seconds\n"
              << "last level: " << LastLevel << '\n';
    mg_.SizeInfo( std::cout);
}

} // end of namespace DROPS

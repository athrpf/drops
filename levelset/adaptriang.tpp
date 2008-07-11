//**************************************************************************
// File:    adaptriang.tpp                                                 *
// Content: adaptive triangulation based on position of the interface      *
//          provided by the levelset function                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

namespace DROPS
{

template <class DistFctT>
void AdapTriangCL::MakeInitialTriang (DistFctT& Dist)
{
    TimerCL time;

    const Uint min_ref_num= f_level_ - c_level_;
    Uint i;
    for (i=0; i<2*min_ref_num; ++i)
        ModifyGridStep( Dist);
    time.Stop();
    std::cout << "MakeInitialTriang: " << i
              << " refinements in " << time.GetTime() << " seconds\n"
              << "last level: " << mg_.GetLastLevel() << '\n';
    mg_.SizeInfo( std::cout);
}

template <class DistFctT>
bool AdapTriangCL::ModifyGridStep (DistFctT& Dist)
// One step of grid change; returns true if modifications were necessary,
// false, if nothing changed.
{
    bool modified= false;
    for (MultiGridCL::TriangTetraIteratorCL it= mg_.GetTriangTetraBegin(),
         end= mg_.GetTriangTetraEnd(); it!=end; ++it)
    {
        double d= 1e99;
        int num_pos= 0;
        for (Uint j=0; j<4; ++j)
        {
            const double dist= GetValue( Dist, *it->GetVertex( j));
            if (dist>=0) ++num_pos;
            d= std::min( d, std::abs( dist));
        }
        for (Uint j=0; j<6; ++j)
        {
            const double dist= GetValue( Dist, *it->GetEdge( j));
            if (dist>=0) ++num_pos;
            d= std::min( d, std::abs( dist));
        }
        d= std::min( d, std::abs( GetValue( Dist, *it)));

        const bool vzw= num_pos!=0 && num_pos!=10; // change of sign
        const Uint l= it->GetLevel();
        // In the shell:      level should be f_level_.
        // Outside the shell: level should be c_level_.
        const Uint soll_level= (d<=width_ || vzw) ? f_level_ : c_level_;

        if (l !=  soll_level || (l == soll_level && !it->IsRegular()) )
        { // tetra will be marked for refinement/removement
            modified= true;
            if (l <= soll_level)
                it->SetRegRefMark();
            else // l > soll_level
                it->SetRemoveMark();
        }
    }
    if (modified) {
        notify_pre_refine();
        mg_.Refine();
        notify_post_refine();
    }
    return modified;
}

inline
void AdapTriangCL::UpdateTriang (const LevelsetP2CL& lset)
{
    TimerCL time;
    modified_= false;
    const int min_ref_num= f_level_ - c_level_;
    int i;
    LevelsetP2CL::const_DiscSolCL sol( lset.GetSolution());

    for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
        (*obs)->pre_refine_sequence();
    for (i= 0; i < 2*min_ref_num; ++i) {
        if (!ModifyGridStep( sol))
            break;
        modified_= true;
    }
    for (ObserverContT::iterator obs= observer_.begin(); obs != observer_.end(); ++obs)
        (*obs)->post_refine_sequence();

    time.Stop();
    std::cout << "UpdateTriang: " << i
              << " refinements/interpolations in " << time.GetTime() << " seconds\n"
              << "last level: " << mg_.GetLastLevel() << '\n';
    mg_.SizeInfo( std::cout);
}

} // end of namespace DROPS

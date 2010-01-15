/// \file adaptriang.tpp
/// \brief adaptive triangulation based on position of the interface provided by the levelset function
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

namespace DROPS
{

template <class DistFctT>
void AdapTriangCL::MakeInitialTriang( DistFctT& Dist)
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif

    time.Reset();
    time.Start();
    const Uint min_ref_num= f_level_;
    Uint i;
    bool modified= true;
    for (i=0; i<2*min_ref_num && modified; ++i)
        modified=ModifyGridStep( Dist);

    time.Stop();
    const double duration=time.GetTime();

    IF_MASTER{
        std::cout << "MakeInitialTriang: " << i
                  << " refinements in " << duration << " seconds\n"
                  << "last level: " << mg_.GetLastLevel() << '\n';
#ifdef _PAR
        DROPS_LOGGER_SETVALUE("MakeInitialTriang",duration);
#endif
    }
    mg_.SizeInfo( std::cout);
}

#ifndef _PAR
template <class DistFctT>
  bool AdapTriangCL::ModifyGridStep( DistFctT& Dist)
#else
template <class DistFctT>
  bool AdapTriangCL::ModifyGridStep( DistFctT& Dist, bool lb)
#endif
/** One step of grid change; returns true if modifications were necessary,
    false, if nothing changed. */
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
#ifdef _PAR
    modified= ProcCL::GlobalOr(modified);
#endif
    if (modified) {
        notify_pre_refine();
        mg_.Refine();
#ifdef _PAR
        pmg_->HandleUnknownsAfterRefine();
        if (lb)
            lb_.DoMigration();
#endif
        notify_post_refine();
#ifdef _PAR
        Assert(!DDD_ConsCheck(), DROPSErrCL("AdapTriangCL::ModifyGridStep: Failure in DDD_ConsCheck"),
               DebugParallelC|DebugParallelNumC|DebugLoadBalC);
#endif
    }
    return modified;
}

inline
void AdapTriangCL::UpdateTriang (const LevelsetP2CL& lset)
/** This function updates the triangulation according to the position of the
    interface provided by the levelset function. Therefore this function marks
    and refines tetras and balance the number of tetras over the processors.
    Also the numerical datas are interpolated to the new triangulation. */
{
#ifndef _PAR
    TimerCL time;
#else
    ParTimerCL time;
#endif
    double duration;

    modified_= false;
    const int min_ref_num= f_level_ - c_level_;
    int i;
    LevelsetP2CL::const_DiscSolCL sol( lset.GetSolution());

    notify_pre_refine_sequence();
    for (i= 0; i < 2*min_ref_num; ++i) {
#ifndef _PAR
        if (!ModifyGridStep( sol))
            break;
#else
        if (!ModifyGridStep(sol, i==2*min_ref_num-1))
            break;
#endif
        modified_= true;
    }
    notify_post_refine_sequence();

    time.Stop();
    duration= time.GetTime();
    IF_MASTER{
        std::cout << "UpdateTriang: " << i
                  << " refinements/interpolations in " << duration << " seconds\n"
                  << "last level: " << mg_.GetLastLevel() << '\n';
#ifdef _PAR
        DROPS_LOGGER_SETVALUE("UpdateTriang",duration);
#endif
    }
    mg_.SizeInfo( std::cout);
}

} // end of namespace DROPS

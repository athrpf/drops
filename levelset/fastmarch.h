/// \file fastmarch.h
/// \brief fast marching method for reparametrization
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


#ifndef DROPS_FASTMARCH_H
#define DROPS_FASTMARCH_H

#include "geom/multigrid.h"
#include "misc/problem.h"
#include "num/spmat.h"
#include "num/interfacePatch.h"
#ifdef _PAR
#  include "parallel/interface.h"
#  include "parallel/exchange.h"
#  include "parallel/partime.h"
#endif
#include <set>

namespace DROPS
{
/// \brief Class for performing the fast marching algorithm for reparametrization of the levelset function
/// \todo (merge) It would be much more "beautifull" if v_ is declared as a pointer.
class FastMarchCL
{
  public:
    enum { Finished= 1, Close= 2, Far= 0};
#ifdef _PAR
    /// \brief structure for sending mark and value
    struct CoupMarkValST
    {
        byte mark;
        double val;
    };
#endif
    typedef SArrayCL<IdxT,4>        ReprTetraT;
    typedef std::vector<ReprTetraT> VertexNeighT;

  private:
    MultiGridCL&               MG_;             // multigrid
#ifndef _PAR
    VecDescCL*                 v_;
    const IdxT                 size_;           // number of DoF
    VectorBaseCL<byte>         Typ_;
#else
    static VecDescCL*          v_;              // values of the levelset function
    static VectorCL            tmpv_;           // value while sending types process (static for DDD)
    static VectorBaseCL<byte>  Typ_;
    static VectorBaseCL<byte>  tmpTyp_;         // typ while sending types process (static for DDD)
    static IdxT                offset_;         // my offset for global numbering
    static VectorBaseCL<IdxT>  allExclusive_;   // offset of all processors
    static std::vector<IdxT>   globNumb_;       // mapping local DoF on global DoF number
    static std::vector<IdxT>   locNumb_;        // mapping global DoF number on local exclusive DoF number

    std::vector<bool>          exclusive_;      // Check if an index is exlusive
    IdxT                       globsize_;       // number of all DoF

    VectorBaseCL<VertexNeighT> GlobNeigh_;                          // neighboors for all DoF
    VectorBaseCL<byte>         GlobalTyp_;                          // typs for all DoF
    VectorBaseCL<Point3DCL>    GlobalCoord_;                        // coordinates for all DoF
    VectorBaseCL<IdxT>         allOffset_;                          // offset of procs
    VectorCL                   GlobalV_;                            // all values of the level set function
    VertexNeighT               TetraList_;                          // all tetras on a proc
    ExchangeCL                *ex_;                                 // for accumulation of values

    const IdxT                 size_;                               // number of DoF
#endif
    VectorCL                   CoordZeroLevel_;                     // Coordinates of zero level of levelset function over all process
    VectorCL                   ValueZeroLevel_;                     // Values of zero level of levelset function over all procs

    std::set<IdxT>             Close_;                              // DoF which have a Finished marked DoF as neighboor
    VectorBaseCL<VertexNeighT> neigh_;                              // neighboor tetras of DoF
    VectorBaseCL<Point3DCL>    Coord_;                              // coordinates of DoF
    VectorCL                   Old_;                                // old values
    VectorBaseCL<IdxT>         map_;
    InterfacePatchCL           patch_;

    void   InitClose();                                 // Mark verts next to finished verts as close and update them. Initialize also neigh_
    IdxT   FindTrial() const;                           // Find vert with minimal value in Close_
    void   Update( const IdxT);                         // Update value on a vert and put this vert into Close_
    double CompValueProj( IdxT, int num, const IdxT upd[3]) const;
    inline void Normalize( double& b) const;                 // Normalize b onto unit interval [0,1]
    inline void Normalize( double& b1, double& b2) const;    // Normalize (b1,b2) onto unit triangle
    /// \brief Store all vertices and values around the zero level of the levelset function
    void InitZeroSet(VectorCL&, VectorCL&);
    /// \brief Approximate Distance to zero level
    double MinDist(IdxT);

    // variants for periodic boundaries
    void   InitClosePer();
    IdxT   FindTrialPer() const;
    void   UpdatePer( const IdxT);
    double CompValueProjPer( IdxT, int num, const IdxT upd[3]) const;
    inline IdxT Map( IdxT i) const { return i<size_ ? i: map_[i-size_]; }

  public:
    FastMarchCL( MultiGridCL& mg, VecDescCL& v)
#ifndef _PAR
      : MG_(mg), v_(&v), size_(v.RowIdx->NumUnknowns()), Typ_(Far, size_) {}
#else
      : MG_(mg), ex_( &v.RowIdx->GetEx()), size_(v.RowIdx->NumUnknowns()), Coord_(size_)
    {
        v_= &v;
        Typ_.resize(Far, size_);
    }
#endif

    /// \brief Finds the zero level of the levelset function and handle verts around it.
    void InitZero( bool ModifyZero= true, int method= 0);
    /// \brief Restore signs of the levelset function according to old signs.
    void RestoreSigns();
    /// \brief Reparametrization of the levelset function with the fast marching algorithm.
    /// This function also calls InitZero and RestoreSigns.
    void Reparam( bool ModifyZero= true, int method= 0);

    /// \brief Reparametrize with Euclidian distance to vertices around the zero level
    void ReparamEuklid( bool ModifyZero=true);

    /// \name variants for periodic boundaries
    //@{
    void InitZeroPer( const BndDataCL<>&, bool ModifyZero= true, int method= 0);
    void ReparamPer( const BndDataCL<>&, bool ModifyZero= true);
    //@}

#ifdef _PAR
    void CreateGlobNumb();                                      ///< Create global numbering for fast marching algorithm
    void DistributeFinished();                                  ///< Communicate Finished marked DoF
    void InitLocNeigh();                                        ///< Compute local neighborhood
    void Collect();                                             ///< Send all information to master-process
    void Distribute();                                          ///< Send all information from master-process to all other procs
    void CleanUp();                                             ///< Free not used memory

    // handlers for DDD (so they are static)
    template<class SimplexT>
      static int HandlerFinishedGather(DDD_OBJ, void*);         ///< Gather finished DoF
    template<class SimplexT>
      static int HandlerFinishedScatter(DDD_OBJ, void*);        ///< Scatter finished DoF
    template<class SimplexT>
      static int HandlerGlobDOFGather(DDD_OBJ, void*);          ///< Gather global number of DoFs
    template<class SimplexT>
      static int HandlerGlobDOFScatter(DDD_OBJ, void*);         ///< Scatter global number of DoFs

    inline IdxT GetLocNum(IdxT globNum);                        ///< Get local number of a DoF from global number
    inline IdxT GetGlobNum(IdxT locNum);                        ///< Get global number of a local DoF
    inline bool IsExclusive(IdxT i);                            ///< Check if a DoF is exclusive

    void DistributeZeroLevel();                                 ///< Distribute all DoF marked as Finished to all procs
#endif
};

inline void FastMarchCL::Normalize( double& b) const
// Normalize b onto unit interval [0,1]
{
    if (b<0) b=0;
    else if (b>1) b=1;
}

inline void FastMarchCL::Normalize( double& b1, double& b2) const
// Normalize (b1,b2) onto unit triangle
{
    Normalize(b1); Normalize(b2); // now (b1,b2) is in unit square
    const double a=(b1+b2-1)/2.;
    if (a>0) { b1-=a; b2-=a; } // projection on diagonal y = -x
}


#ifdef _PAR
// Declaration of wrapper for gathering
// and scattering data for giving these to DDD
//--------------------------------------------
extern "C" int HandlerFinishedGatherVertexC(DDD_OBJ objp, void* buf);
extern "C" int HandlerFinishedGatherEdgeC(DDD_OBJ objp, void* buf);
extern "C" int HandlerFinishedScatterVertexC(DDD_OBJ objp, void* buf);
extern "C" int HandlerFinishedScatterEdgeC(DDD_OBJ objp, void* buf);

extern "C" int HandlerGlobDOFGatherVertexC(DDD_OBJ objp, void* buf);
extern "C" int HandlerGlobDOFGatherEdgeC(DDD_OBJ objp, void* buf);
extern "C" int HandlerGlobDOFScatterVertexC(DDD_OBJ objp, void* buf);
extern "C" int HandlerGlobDOFScatterEdgeC(DDD_OBJ objp, void* buf);

// Definition of inline functions
//--------------------------------------------
IdxT FastMarchCL::GetGlobNum(IdxT locNum)
{
    return globNumb_[locNum];
}

bool FastMarchCL::IsExclusive(IdxT i)
{
    return exclusive_[i];
}
#endif  // end of parallel
} // end of namespace DROPS

#endif

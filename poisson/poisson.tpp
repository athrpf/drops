/// \file poisson.tpp
/// \brief classes that constitute the poisson-problem
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Marcus Soemers, Volker Reichelt, Liang Zhang; SC RWTH Aachen: Oliver Fortmeier

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

#include "num/discretize.h"
#include "num/accumulator.h"
#include "poisson/poissonaccus.tpp"
namespace DROPS
{

//========================================================================================================
//
//                       Set up matrices and right hand side (rhs) in poisson P1 problem
//
//========================================================================================================


/// - Setup-System routines


/// \brief This function is called by different setup routines in poissonP1 class
/// Depending on corresponding pointers are NULL or not, this function will decide 
/// which part of system should be assembled in accumulation way
//========================================================================================================
// Explations of part of arguments of SetupPartialSysterm_P1 function
// MatrixCL* Amat       Stiffness matrix  
// MatrixCL* Mmat       mass matrix
// MatrixCL* Umat       convection matrix
// VecDescCL* cplA      boundary coupling term corresponding to Amat
// VecDescCL* cplM      boundary coupling term corresponding to Mmat
// VecDescCL* cplU      boundary coupling term corresponding to Umat
// VecDescCL* f         Source functions
// SUPGCL& supg         SUPG stabilization 
// bool adjoint         Adjoint problem flag        
//========================================================================================================
template<class Coeff>
void SetupPartialSystem_P1( const MultiGridCL& MG, const Coeff& , MatrixCL* Amat, MatrixCL* Mmat, 
                          MatrixCL* Umat, VecDescCL* cplA, VecDescCL* cplM, VecDescCL* cplU, VecDescCL* f, 
                          const BndDataCL<> * BndData_,
                          IdxDescCL& RowIdx, IdxDescCL& ColIdx, double t, SUPGCL& supg, bool adjoint)
{
    //Create tuple
    TetraAccumulatorTupleCL accus;
    //Accumulators
    StiffnessAccumulator_P1CL<Coeff,Quad5CL> * accua = 0;
    MassAccumulator_P1CL<Coeff,Quad5CL> * accum = 0;
    ConvectionAccumulator_P1CL<Coeff,Quad3CL> * accuc = 0;
    SourceAccumulator_P1CL<Coeff,Quad5CL> * accuf = 0;
    
    //register accumulator
    if (Amat !=0 || cplA !=0){
        accua = new StiffnessAccumulator_P1CL<Coeff,Quad5CL> (MG,BndData_,Amat,cplA,RowIdx,ColIdx,supg,t);
        accus.push_back( accua);
    }
    if (Mmat !=0 || cplM !=0){
        accum = new MassAccumulator_P1CL<Coeff,Quad5CL> (MG,BndData_,Mmat,cplM,RowIdx,ColIdx,supg,t);
        accus.push_back( accum);
    }
    if (Umat !=0 || cplU !=0){
        accuc = new ConvectionAccumulator_P1CL<Coeff,Quad3CL> (MG,BndData_,Umat,cplU,RowIdx,ColIdx,t,adjoint);
        accus.push_back( accuc);
    }
    if (f !=0){
        accuf = new SourceAccumulator_P1CL<Coeff,Quad5CL> (MG,BndData_,f,RowIdx,supg,t);
        accus.push_back( accuf);
    }
    
    //run accumulation
    accumulate( accus, MG, RowIdx.TriangLevel(),RowIdx.GetMatchingFunction(), RowIdx.GetBndInfo());  
    if (accua != 0) delete accua;
    if (accum != 0) delete accum;
    if (accuc != 0) delete accuc;
    if (accuf != 0) delete accuf;
}


template<class Coeff>
void PoissonP1CL<Coeff>::SetupSystem(MLMatDescCL& matA, VecDescCL& b) const
///Go throught every level to Setup system in P1
{
    MLMatrixCL::iterator  itA    = matA.Data.begin();
    MLIdxDescCL::iterator itRow  = matA.RowIdx->begin();
    MLIdxDescCL::iterator itCol  = matA.ColIdx->begin();
    for ( size_t lvl=0; lvl < matA.Data.size(); ++lvl, ++itRow, ++itCol, ++itA){
        VecDescCL * rhs = (lvl == matA.Data.size()-1) ? &b : 0;
        SetupPartialSystem_P1(MG_,Coeff_, &*itA,0,0,rhs,0,0,rhs, &BndData_,*itRow,*itCol, /* time */ 0.,supg_, /* adjoint */ false);
    }
}

template<class Coeff>
void PoissonP1CL<Coeff>::SetupInstatSystem( MLMatDescCL& matA, MLMatDescCL& matM, double t) const
///Go throught every multigrid level to Setup left hand side of instationary system in P1
{
    MLIdxDescCL::iterator itRow  = matA.RowIdx->begin();
    MLIdxDescCL::iterator itCol  = matA.ColIdx->begin();
    MLMatrixCL::iterator  itM    = matM.Data.begin();
    for ( MLMatrixCL::iterator itA= matA.Data.begin(); itA != matA.Data.end(); ++itA, ++itM, ++itRow, ++itCol)
        SetupPartialSystem_P1(MG_,Coeff_,&*itA, &*itM,0,0,0,0,0,0,*itRow, *itCol,t,supg_,false);
}

template<class Coeff>
void PoissonP1CL<Coeff>::SetupConvection( MLMatDescCL& matU, VecDescCL& vU, double t) const
///Go throught every multigrid level to setup convection 
{
    MLMatrixCL::iterator  itU    = matU.Data.begin();
    MLIdxDescCL::iterator itRow  = matU.RowIdx->begin();
    MLIdxDescCL::iterator itCol  = matU.ColIdx->begin();
    for ( size_t lvl=0; lvl < matU.Data.size(); ++lvl, ++itRow, ++itCol, ++itU)
        SetupPartialSystem_P1(MG_,Coeff_,0, 0,&*itU,0,0,(lvl == matU.Data.size()-1) ? &vU : 0,0,&BndData_,*itRow, *itCol,t,supg_,adjoint_);
        
}

template<class Coeff>
void PoissonP1CL<Coeff>::SetupInstatRhs(VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double) const
///Setup right hand side of instationary system in P1 only last level
{
  SetupPartialSystem_P1(MG_,Coeff_,0,0,0,&vA,&vM,0,&vf,&BndData_,*vA.RowIdx,*vA.RowIdx,tA,supg_,false);
}

template<class Coeff>
void PoissonP1CL<Coeff>::SetNumLvl( size_t n)
{
    match_fun match= MG_.GetBnd().GetMatchFun();
    idx.resize( n, P1_FE, BndData_, match);
    A.Data.resize( idx.size());
    M.Data.resize( idx.size());
    U.Data.resize( idx.size());
}

template <class Coeff>
void PoissonP1CL<Coeff>::Init( VecDescCL& vec, instat_scalar_fun_ptr func, double t0) const
///Setup initial condition for instationary problem
{
    Uint lvl= vec.GetLevel(),
         idx= vec.RowIdx->GetIdx();


    for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send= const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (sit->Unknowns.Exist(idx))
        {
            vec.Data[sit->Unknowns(idx)]= func( sit->GetCoord(), t0);
        }
    }

}



template <class Coeff>
void PoissonP1CL<Coeff>::SetupGradSrc(VecDescCL& src, instat_scalar_fun_ptr T, instat_scalar_fun_ptr dalpha, double t) const
///Special rhs for IA2 sensitivity problem
{
  src.Clear( t);
  const Uint lvl = src.GetLevel(),
             idx = src.RowIdx->GetIdx();
  Point3DCL G[4];

  double det;
  double absdet;
  IdxT UnknownIdx[4];
  Quad2CL<> rhs, quad_a;

  for (MultiGridCL::const_TriangTetraIteratorCL
    sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl),
    send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
    sit != send; ++sit)
  {
    P1DiscCL::GetGradients(G,det,*sit);
    absdet= std::fabs(det);

    quad_a.assign( *sit, dalpha, t);
    const double int_a= quad_a.quad( absdet);
    Point3DCL gradT;

    for(int i=0; i<4; ++i)
    {
      gradT+= G[i]*T(sit->GetVertex(i)->GetCoord(), t);
      UnknownIdx[i]= sit->GetVertex(i)->Unknowns.Exist(idx) ? sit->GetVertex(i)->Unknowns(idx)
                                                            : NoIdx;
    }

    for(int i=0; i<4;++i)    // assemble row i
    {
      if (sit->GetVertex(i)->Unknowns.Exist(idx)) // vertex i is not on a Dirichlet boundary
      {
        src.Data[UnknownIdx[i]]-= int_a*inner_prod( gradT, G[i]);
        if ( BndData_.IsOnNatBnd(*sit->GetVertex(i)) )
          for (int f=0; f < 3; ++f)
            if ( sit->IsBndSeg(FaceOfVert(i, f)) )
            {
              Point3DCL n;
              sit->GetOuterNormal(FaceOfVert(i, f), n);
              src.Data[UnknownIdx[i]]+=
                P1DiscCL::Quad2D(*sit, FaceOfVert(i, f), dalpha, i,  t) * inner_prod( gradT, n);
            }
      }
    }
  }
}


//=======================================================================================================
//
//            check computed solution, estimate error etc. in poisson P1 problem
//
//========================================================================================================
template <class Coeff>
double PoissonP1CL<Coeff>::CheckSolution(const VecDescCL& lsg,
  instat_scalar_fun_ptr Lsg, double t) const
{
  double diff, maxdiff=0, norm2= 0, L2=0;
  Uint lvl=lsg.GetLevel(),
       Idx=lsg.RowIdx->GetIdx();

  const_DiscSolCL sol(&lsg, &GetBndData(), &GetMG());

  std::cout << "Difference to exact solution:" << std::endl;

  for (MultiGridCL::const_TriangTetraIteratorCL
    sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl),
    send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
    sit != send; ++sit)
  {
    double absdet= sit->GetVolume()*6.,
           sum= 0;

    for(Uint i=0; i<4; ++i)
    {
      diff= (sol.val(*sit->GetVertex(i)) - Lsg(sit->GetVertex(i)->GetCoord(),t));
      sum+= diff*diff;
    }
    sum/= 120;
    diff= sol.val(*sit, 0.25, 0.25, 0.25) - Lsg(GetBaryCenter(*sit),t);
    sum+= 2./15. * diff*diff;
    L2+= sum*absdet;
  }
#ifdef _PAR
  L2= ProcCL::GlobalSum(L2);
#endif
  L2= std::sqrt(L2);

  for (MultiGridCL::const_TriangVertexIteratorCL
    sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl),
    send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
    sit != send; ++sit)
  {
    if (sit->Unknowns.Exist(Idx))
#ifdef _PAR
      if (sit->IsExclusive())
#endif
      {
        diff= std::fabs( Lsg(sit->GetCoord(),t) - lsg.Data[sit->Unknowns(Idx)] );
        norm2+= diff*diff;
        if (diff>maxdiff)
        {
          maxdiff= diff;
        }
      }
  }

  int Lsize = lsg.Data.size();
#ifdef _PAR
  Lsize   = lsg.RowIdx->GetGlobalNumUnknowns(MG_);
  norm2   = ProcCL::GlobalSum(norm2, Drops_MasterC);
  maxdiff = ProcCL::GlobalMax(maxdiff, Drops_MasterC);
#endif
  std::cout << "  2-Norm= " << std::sqrt(norm2)
            << "\nw-2-Norm= " << std::sqrt(norm2/Lsize)
            << "\nmax-Norm= " << maxdiff
            << "\n L2-Norm= " << L2 << std::endl;
  return L2;
}

//=======================================================================================================
//
//            error estimation in poisson P1 problem
//
//========================================================================================================

template<class Coeff>
void PoissonP1CL<Coeff>::GetDiscError(const MLMatDescCL& A, instat_scalar_fun_ptr Lsg) const
{
    Uint lvl= A.GetColLevel(),
         idx= A.ColIdx->GetIdx();
    VectorCL lsg(A.Data.num_cols());

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (sit->Unknowns.Exist(idx))
        {
            lsg[sit->Unknowns(idx)]= Lsg(sit->GetCoord(), 0.0);
        }
    }

    VectorCL res( A.Data*lsg-b.Data);
    std::cout <<"|| Ax - b || = "<< norm( res)<<", max "<< supnorm( res)
              <<" (with x=continuous solution)"<<std::endl;
}



template<class Coeff>
bool PoissonP1CL<Coeff>::EstimateError (const VecDescCL& lsg,
   const double rel_loc_tol, double& globalerr, est_fun est)
{
    const Uint lvl= lsg.GetLevel();
    Uint num_ref= 0;
    Uint num_un_ref= 0;
    const Uint num_tetra= std::distance(MG_.GetTriangTetraBegin(lvl), MG_.GetTriangTetraEnd(lvl));
    const double exp_err= globalerr/0.875; // divident is the volume of the domain
    double localerr;
    double localerr_dist;

    globalerr= 0.0;
    for (MultiGridCL::TriangTetraIteratorCL sit=MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        localerr= est(*sit, lsg, GetBndData());
        localerr_dist= localerr/std::sqrt(sit->GetVolume());
        globalerr+= localerr*localerr;
        if ( sit->IsUnrefined() ) {
            if (localerr_dist>rel_loc_tol*exp_err)
            {
                sit->SetRegRefMark();
                ++num_ref;
            }
            else if (localerr_dist<exp_err/rel_loc_tol)
                {
                    sit->SetRemoveMark();
                    ++num_un_ref;
                }
        }
    }
    globalerr= std::sqrt(globalerr);
    std::cout << "Estimated global W1_2-error: " << globalerr << ". Marked " << num_ref << ", unmarked " << num_un_ref
              << " out of " << num_tetra << " tetrahedrons." << std::endl;
    return num_ref || num_un_ref;
}


template<class Coeff>
double PoissonP1CL<Coeff>::ResidualErrEstimator(const TetraCL& t, const VecDescCL& sol, const BndDataCL& Bnd)
// Based on R. Verfuerth, "A review of a posteriori error estimation and adaptive
// mesh refinement techniques"
// chapter 1.2
{
    double _err= 0.0;
    Point3DCL cc_center;  // center of circum-circle
    double cc_radius;     // radius of circum-circle
    double cc_rad_Face;
    Uint trilevel= sol.GetLevel();
    const Uint Idx= sol.RowIdx->GetIdx();

    double det;
    SMatrixCL<3,4> H;
    // Gradient of hat-function for vertex i is in col i
    P1DiscCL::GetGradients(H,det,t);
    const double absdet= std::fabs(det);
    // Integral over tetra
    circumcircle(t, cc_center, cc_radius);
    const double P0f= Quad3PosWeightsCL::Quad(t, CoeffCL::f, 0.0)*6.; // absdet/vol = 6.
    _err+= 4.*cc_radius*cc_radius*P0f*P0f*absdet/6.;

    // Integrals over boundary - distinguish between inner-boundary and domain-boundary
    Point3DCL normal; // normal of face
    double dir;       // 1.0, if normal points outward, -1.0 otherwise
    for (Uint face=0; face<NumFacesC; ++face)
    {
        circumcircle(t, face, cc_center, cc_rad_Face);
        if ( t.IsBndSeg(face) )
        {
            t.GetOuterNormal(face, normal);
            const BndIdxT bidx= t.GetBndIdx(face);
            const typename BndDataCL::BndSegT bdat= Bnd.GetBndSeg( bidx);
            if ( bdat.IsNatural() )
            {
                Point3DCL vc3D[3];
                const VertexCL* v[4];

                double n_du= 0.0;  // scalar product of Du and the normal
                for (int i=0; i<4; ++i)
                {
                    v[i]= t.GetVertex(i);
                    n_du+= (v[i]->Unknowns.Exist(Idx) ? sol.Data[v[i]->Unknowns(Idx)]
                               : Bnd.GetDirBndValue(*v[i]))
                          *(H(0,i)*normal[0] + H(1,i)*normal[1] + H(2,i)*normal[2]);
                }
                for (int i=0; i<3; ++i)
                    vc3D[i]= v[VertOfFace(face,i)]->GetCoord();

                const double f0= bdat.GetBndVal(vc3D[0]) - n_du;
                const double f1= bdat.GetBndVal(vc3D[1]) - n_du;
                const double f2= bdat.GetBndVal(vc3D[2]) - n_du;
                const double fb= bdat.GetBndVal( 1./3.*(vc3D[0] + vc3D[1] + vc3D[2]) ) - n_du;    //Barycenter of Face
                const double absdet= FuncDet2D(v[VertOfFace(face,1)]->GetCoord() - v[0]->GetCoord(), v[2]->GetCoord() - v[0]->GetCoord());
                _err+= 2.*cc_rad_Face * ( (f0*f0 + f1*f1 + f2*f2)/24. + 3./8.*fb*fb )*absdet/2.;
            }
            // TODO: How do we handle non-null Dirichlet-boundary-data
        }
        else
        {
            const double absdet2D= t.GetNormal(face, normal, dir);
            const TetraCL& neigh= *t.GetNeighInTriang(face, trilevel);
            double ndet;
            SMatrixCL<3,4> nH;
            P1DiscCL::GetGradients(nH, ndet, neigh);
            const VertexCL* v[4];
            const VertexCL* nv;
            double jump= 0.0;
            for (int i=0; i<4; ++i)
            {
                v[i]= t.GetVertex(i);
                nv= neigh.GetVertex(i);
                jump-= dir
                      *(v[i]->Unknowns.Exist(Idx) ? sol.Data[v[i]->Unknowns(Idx)]
                        : Bnd.GetDirBndValue(*v[i]) )
                      *(H(0,i)*normal[0] + H(1,i)*normal[1] + H(2,i)*normal[2]);
                jump+= dir
                      *(nv->Unknowns.Exist(Idx) ? sol.Data[nv->Unknowns(Idx)]
                        : Bnd.GetDirBndValue(*nv) )
                      *(nH(0,i)*normal[0] + nH(1,i)*normal[1] + nH(2,i)*normal[2]);
            }
            _err+= /* 0.5 * 2.* */cc_rad_Face * jump*jump * absdet2D/2.;
        }
    }
    return _err;
}


template<class Coeff>
double PoissonP1CL<Coeff>::ResidualErrEstimatorL2(const TetraCL& t, const VecDescCL& sol, const BndDataCL& Bnd)
// Based on R. Verfuerth, "A review of a posteriori error estimation and adaptive
// mesh refinement techniques"
// chapter 3.2
{
    double _err= 0.0;
    Point3DCL cc_center;  // center of circum-circle
    double cc_radius;     // radius of circum-circle
    double cc_rad_Face;
    Uint trilevel= sol.GetLevel();
    const Uint Idx= sol.RowIdx->GetIdx();

    double det;
    SMatrixCL<3,4> H;
    // Gradient of hat-function for vertex i is in col i
    P1DiscCL::GetGradients(H,det,t);
    const double absdet= std::fabs(det);
    // Integral over tetra
    circumcircle(t, cc_center, cc_radius);
    const double P0f= Quad3PosWeightsCL::Quad(t, &CoeffCL::f, 0.0)*6.; // absdet/vol = 6.
    _err+= 4.*cc_radius*cc_radius*P0f*P0f*absdet/6.;

    // Integrals over boundary - distinguish between inner-boundary and domain-boundary
    Point3DCL normal; // normal of face
    double dir;       // 1.0, if normal points outward, -1.0 otherwise
    for (Uint face=0; face<NumFacesC; ++face)
    {
        circumcircle(t, face, cc_center, cc_rad_Face);
        if ( t.IsBndSeg(face) )
        {
            t.GetOuterNormal(face, normal);
            const BndIdxT bidx= t.GetBndIdx(face);
            const typename BndDataCL::BndSegT bdat= Bnd.GetBndSeg( bidx);
            if ( bdat.IsNatural() )
            {
                Point3DCL vc3D[3];
                const VertexCL* v[4];

                double n_du= 0.0;  // scalar product of Du and the normal
                for (int i=0; i<4; ++i)
                {
                    v[i]= t.GetVertex(i);
                    n_du+= (v[i]->Unknowns.Exist(Idx) ? sol.Data[v[i]->Unknowns(Idx)]
                               : Bnd.GetDirBndValue(*v[i]) )
                          *(H(0,i)*normal[0] + H(1,i)*normal[1] + H(2,i)*normal[2]);
                }
                for (int i=0; i<3; ++i)
                    vc3D[i]= v[VertOfFace(face,i)]->GetCoord();

                const double f0= bdat.GetBndVal(vc3D[0]) - n_du;
                const double f1= bdat.GetBndVal(vc3D[1]) - n_du;
                const double f2= bdat.GetBndVal(vc3D[2]) - n_du;
                const double fb= bdat.GetBndVal( 1./3.*(vc3D[0] + vc3D[1] + vc3D[2]) ) - n_du;    //Barycenter of Face
                const double absdet= FuncDet2D(v[VertOfFace(face,1)]->GetCoord() - v[0]->GetCoord(), v[2]->GetCoord() - v[0]->GetCoord());
                _err+= 2.*cc_rad_Face * ( (f0*f0 + f1*f1 + f2*f2)/24. + 3./8.*fb*fb )*absdet/2.;
            }
            // TODO: How do we handle non-null Dirichlet-boundary-data
        }
        else
        {
            const double absdet2D= t.GetNormal(face, normal, dir);
            const TetraCL& neigh= *t.GetNeighInTriang(face, trilevel);
            double ndet;
            SMatrixCL<3,4> nH;
            P1DiscCL::GetGradients(nH, ndet, neigh);
            const VertexCL* v[4];
            const VertexCL* nv;
            double jump= 0.0;
            for (int i=0; i<4; ++i)
            {
                v[i]= t.GetVertex(i);
                nv= neigh.GetVertex(i);
                jump-= dir
                      *(v[i]->Unknowns.Exist(Idx) ? sol.Data[v[i]->Unknowns(Idx)]
                        : Bnd.GetDirBndValue(*v[i]) )
                      *(H(0,i)*normal[0] + H(1,i)*normal[1] + H(2,i)*normal[2]);
                jump+= dir
                      *(nv->Unknowns.Exist(Idx) ? sol.Data[nv->Unknowns(Idx)]
                        : Bnd.GetDirBndValue(*nv) )
                      *(nH(0,i)*normal[0] + nH(1,i)*normal[1] + nH(2,i)*normal[2]);
            }
            _err+= /* 0.5 * 2.* */cc_rad_Face * jump*jump * absdet2D/2.;
        }
    }
    return 4.*cc_radius*cc_radius*_err;
}

//=======================================================================================================
//
//            Local quadrature functions in poisson P2 problem
//            To do: replace them
//
//========================================================================================================

// Copy of functions GetGradientsOnRef, MakeGradients,
// Quad(t,sf,i,j), Quad(t,sf,i)(new!!!) and QuadGrad (for P2P1)
// from stokes-Verzeichnis (stoks.tpp)

inline void GetGradientsOnRef(SMatrixCL<3,5>* GRef)
{
    SVectorCL<3> vec;

    for(int i=0; i<10; ++i)
    {
        vec= FE_P2CL::DHRef(i,0,0,0);                // vertex 0
        GRef[i](0,0)= vec[0];   GRef[i](1,0)= vec[1];   GRef[i](2,0)= vec[2];
        vec= FE_P2CL::DHRef(i,1,0,0);                // vertex 1
        GRef[i](0,1)= vec[0];   GRef[i](1,1)= vec[1];   GRef[i](2,1)= vec[2];
        vec= FE_P2CL::DHRef(i,0,1,0);                // vertex 2
        GRef[i](0,2)= vec[0];   GRef[i](1,2)= vec[1];   GRef[i](2,2)= vec[2];
        vec= FE_P2CL::DHRef(i,0,0,1);                // vertex 3
        GRef[i](0,3)= vec[0];   GRef[i](1,3)= vec[1];   GRef[i](2,3)= vec[2];
        vec= FE_P2CL::DHRef(i,1./4.,1./4.,1./4.);    // barycenter
        GRef[i](0,4)= vec[0];   GRef[i](1,4)= vec[1];   GRef[i](2,4)= vec[2];
    }
}


inline void MakeGradients (SMatrixCL<3,5>* G, const SMatrixCL<3,5>* GRef, const SMatrixCL<3,3>& T)
{
    for(int i=0; i<10; ++i)
        G[i]= T*GRef[i];
}


/// \brief quadratur of order 2 for \f$ \int f\phi_i\phi_j \f$ for quadratic hat functions
inline double Quad( const TetraCL& s, instat_scalar_fun_ptr f, int i, int j, double t= 0.0)
{
    double a[5];
    if (i>j) std::swap(i,j);
    switch(i*10+j)
    {
      case  0: a[0]= 1./1260.; a[1]= a[2]= a[3]= 0.; a[4]= 1./630.; break;
      case  1: a[0]= a[1]= -1./8505.; a[2]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case  2: a[0]= a[2]= -1./8505.; a[1]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case  3: a[0]= a[3]= -1./8505.; a[1]= a[2]= 11./136080.; a[4]= 4./8505; break;
      case  4: a[0]= 1./2520.; a[1]= -1./2520.; a[2]= a[3]= 0.; a[4]= -1./630.; break;
      case  5: a[0]= 1./2520.; a[2]= -1./2520.; a[1]= a[3]= 0.; a[4]= -1./630.; break;
      case  6: a[0]= a[3]= 1./8505.; a[1]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case  7: a[0]= 1./2520.; a[3]= -1./2520.; a[1]= a[2]= 0.; a[4]= -1./630.; break;
      case  8: a[0]= a[2]= 1./8505.; a[1]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case  9: a[0]= a[1]= 1./8505.; a[2]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 11: a[1]= 1./1260.; a[0]= a[2]= a[3]= 0.; a[4]= 1./630.; break;
      case 12: a[1]= a[2]= -1./8505.; a[0]= a[3]= 11./136080.; a[4]= 4./8505; break;
      case 13: a[1]= a[3]= -1./8505.; a[0]= a[2]= 11./136080.; a[4]= 4./8505; break;
      case 14: a[1]= 1./2520.; a[0]= -1./2520.; a[2]= a[3]= 0.; a[4]= -1./630.; break;
      case 15: a[1]= a[3]= 1./8505.; a[0]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 16: a[1]= 1./2520.; a[2]= -1./2520.; a[0]= a[3]= 0.; a[4]= -1./630.; break;
      case 17: a[1]= a[2]= 1./8505.; a[0]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 18: a[1]= 1./2520.; a[3]= -1./2520.; a[0]= a[2]= 0.; a[4]= -1./630.; break;
      case 19: a[0]= a[1]= 1./8505.; a[2]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 22: a[2]= 1./1260.; a[0]= a[1]= a[3]= 0.; a[4]= 1./630.; break;
      case 23: a[2]= a[3]= -1./8505.; a[0]= a[1]= 11./136080.; a[4]= 4./8505; break;
      case 24: a[2]= a[3]= 1./8505.; a[0]= a[1]= -19./68040.; a[4]= -1./486.; break;
      case 25: a[2]= 1./2520.; a[0]= -1./2520.; a[1]= a[3]= 0.; a[4]= -1./630.; break;
      case 26: a[2]= 1./2520.; a[1]= -1./2520.; a[0]= a[3]= 0.; a[4]= -1./630.; break;
      case 27: a[2]= a[1]= 1./8505.; a[0]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 28: a[2]= a[0]= 1./8505.; a[1]= a[3]= -19./68040.; a[4]= -1./486.; break;
      case 29: a[2]= 1./2520.; a[3]= -1./2520.; a[0]= a[1]= 0.; a[4]= -1./630.; break;
      case 33: a[3]= 1./1260.; a[0]= a[1]= a[2]= 0.; a[4]= 1./630.; break;
      case 34: a[3]= a[2]= 1./8505.; a[0]= a[1]= -19./68040.; a[4]= -1./486.; break;
      case 35: a[3]= a[1]= 1./8505.; a[0]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 36: a[3]= a[0]= 1./8505.; a[1]= a[2]= -19./68040.; a[4]= -1./486.; break;
      case 37: a[3]= 1./2520.; a[0]= -1./2520.; a[1]= a[2]= 0.; a[4]= -1./630.; break;
      case 38: a[3]= 1./2520.; a[1]= -1./2520.; a[0]= a[2]= 0.; a[4]= -1./630.; break;
      case 39: a[3]= 1./2520.; a[2]= -1./2520.; a[0]= a[1]= 0.; a[4]= -1./630.; break;
      case 44: a[0]= a[1]= 37./17010.; a[2]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 45: a[0]= 1./972.; a[1]= a[2]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 46: a[1]= 1./972.; a[0]= a[2]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 47: a[0]= 1./972.; a[1]= a[3]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 48: a[1]= 1./972.; a[0]= a[3]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 49: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 55: a[0]= a[2]= 37./17010.; a[1]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 56: a[2]= 1./972.; a[0]= a[1]= 2./8505.; a[3]= -19./34020.; a[4]= 46./8505.; break;
      case 57: a[0]= 1./972.; a[2]= a[3]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 58: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 59: a[2]= 1./972.; a[0]= a[3]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 66: a[1]= a[2]= 37./17010.; a[0]= a[3]= -17./17010.; a[4]= 88./8505.; break;
      case 67: a[0]= a[1]= a[2]= a[3]= 1./11340.; a[4]= 8./2835.; break;
      case 68: a[1]= 1./972.; a[2]= a[3]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 69: a[2]= 1./972.; a[1]= a[3]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 77: a[0]= a[3]= 37./17010.; a[1]= a[2]= -17./17010.; a[4]= 88./8505.; break;
      case 78: a[3]= 1./972.; a[0]= a[1]= 2./8505.; a[2]= -19./34020.; a[4]= 46./8505.; break;
      case 79: a[3]= 1./972.; a[0]= a[2]= 2./8505.; a[1]= -19./34020.; a[4]= 46./8505.; break;
      case 88: a[1]= a[3]= 37./17010.; a[0]= a[2]= -17./17010.; a[4]= 88./8505.; break;
      case 89: a[3]= 1./972.; a[1]= a[2]= 2./8505.; a[0]= -19./34020.; a[4]= 46./8505.; break;
      case 99: a[2]= a[3]= 37./17010.; a[0]= a[1]= -17./17010.; a[4]= 88./8505.; break;
      default: throw DROPSErrCL("Quad(i,j): no such shape function");
    }
    double sum= a[4]*f(GetBaryCenter(s), t);
    for(Uint i=0; i<4; ++i)
        sum+= a[i]*f(s.GetVertex(i)->GetCoord(), t);
    return sum;
}


inline double Quad( const TetraCL& s, instat_scalar_fun_ptr coeff, int i, double t= 0.0)
{
    double f[5];

    if (i<4) // hat function on vert
    {
        f[0]= coeff( s.GetVertex(i)->GetCoord(), t);
        for (int k=0, l=1; k<4; ++k)
            if (k!=i) f[l++]= coeff( s.GetVertex(k)->GetCoord(), t );
        f[4]= coeff( GetBaryCenter(s), t);
        return f[0]/504. - (f[1] + f[2] + f[3])/1260. - f[4]/126.;
    }
    else  // hat function on edge
    {
        const double ve= 4./945.,  // coeff for verts of edge
                     vn= -1./756.,  // coeff for other verts
                     vs= 26./945.;   // coeff for barycenter
        double a[4];
        a[VertOfEdge(i-4,0)]= a[VertOfEdge(i-4,1)]= ve;
        a[VertOfEdge(OppEdge(i-4),0)]= a[VertOfEdge(OppEdge(i-4),1)]= vn;

        double sum= vs * coeff( GetBaryCenter(s), t );
        for(int k=0; k<4; ++k)
            sum+= a[k] * coeff( s.GetVertex(k)->GetCoord(), t );

        return sum;
    }
}

//Neumann boundary condition for P2 problem
inline double QuadP2( const TetraCL& sit, Uint face, const BndDataCL<>& BndData_, int m)
{
        double NeumannBnd;
        LocalP2CL<>p2[10];
        for (int k=0; k<10; ++k)
            p2[k][k]=1.;
        BaryCoordCL BaryC_face[3];
        // std::memset( BaryC_face, 0, 3*sizeof( BaryCoordCL));
        const VertexCL* v[3];
        for (Uint n= 0; n < 3; ++n) 
        v[n]= sit.GetVertex( VertOfFace( face, n));
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<4;j++)
              if(sit.GetVertex(j)==v[i])
                 BaryC_face[i][j]=1.0;
        }
        Quad5_2DCL<> bnd, P2[10];
        P2[m].assign(p2[m], BaryC_face);
        bnd.assign(sit, BaryC_face, BndData_.GetBndFun(sit.GetBndIdx(face)));
        const double absdet_face=FuncDet2D(v[1]->GetCoord()-v[0]->GetCoord(),
                                           v[2]->GetCoord()-v[0]->GetCoord());  
        NeumannBnd=Quad5_2DCL<>(bnd*P2[m]).quad(absdet_face);
        return NeumannBnd;
}


inline double QuadGrad(const SMatrixCL<3,5>* G, int i, int j)
// computes int( grad(phi_i) * grad(phi_j) ) for P2-elements on ref. tetra
{
    SVectorCL<5> tmp(0.0);

    for(int k=0; k<5; ++k)
        for(int l=0; l<3; ++l)
            tmp[k]+= G[i](l,k) * G[j](l,k);

    return ( tmp[0] + tmp[1] + tmp[2] + tmp[3] )/120. + 2./15.*tmp[4];
}


//=======================================================================================================
//
//            Setup matrices and right hand side in poisson P2 problem
//
//========================================================================================================


template<class Coeff>
void SetupSystem_P2( const MultiGridCL& MG, const Coeff&, const BndDataCL<> BndData_, MatrixCL& Amat, VecDescCL* b, IdxDescCL& RowIdx, IdxDescCL& ColIdx)
// Sets up the stiffness matrix and right hand side
{
    if (b!=0) b->Clear( 0.0);

    const Uint lvl    = RowIdx.TriangLevel(),
               idx    = RowIdx.GetIdx();
    MatrixBuilderCL A( &Amat, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());

    IdxT Numb[10];
    bool IsOnDirBnd[10];

#ifndef _PAR
    std::cout << "entering SetupSystem: " << ColIdx.NumUnknowns()
              << " unknowns, " << std::endl;
#endif

// fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup[10][10];
    double det, absdet;
    double tmp;

    GetGradientsOnRef(GradRef);

    for (MultiGridCL::const_TriangTetraIteratorCL sit=MG.GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= std::fabs(det);

        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= BndData_.IsOnDirBnd( *sit->GetVertex(i))))
                Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        }

        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= BndData_.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }

        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // negative dot-product of the gradients
                coup[i][j] = Coeff::alpha * QuadGrad( Grad, i, j)*absdet;
                coup[i][j]+= Quad(*sit, Coeff::q, i, j, 0.0)*absdet;
                coup[j][i] = coup[i][j];
            }

        // assemble

        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        A ( Numb[i], Numb[j] ) += coup[j][i];
                    }
                    else // coupling with vert/edge j on right-hand-side
                        if (b!=0)
                        {
                            tmp= j<4 ? BndData_.GetDirBndValue(*sit->GetVertex(j), 0.0)
                                    : BndData_.GetDirBndValue(*sit->GetEdge(j-4), 0.0);
                            b->Data[Numb[i]]-= coup[j][i] * tmp;
                        }
                }
                if (b!=0)
                {
                    tmp= Quad(*sit, Coeff::f, i, 0.0)*absdet;
                    b->Data[Numb[i]]+= tmp;

                    if ( i<4 ? BndData_.IsOnNatBnd(*sit->GetVertex(i))
                            : BndData_.IsOnNatBnd(*sit->GetEdge(i-4)) ) // vert/edge i is on natural boundary
                    {
                      Uint face;
                      int n=i<4? 3:2;  //Three faces for a vertex and two faces for a point on the edge;
                      for (int f=0; f < n; ++f)
                        {
                         face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);  
                         if ( sit->IsBndSeg(face))
                            {
                             b->Data[Numb[i]]+=QuadP2(*sit, face, BndData_, i);
                            }
                        }
                    }
                }
            }
    }
#ifndef _PAR
    std::cout << "done: value part fill" << std::endl;
#endif

    A.Build();
}

template<class Coeff>
void PoissonP2CL<Coeff>::SetupSystem(MLMatDescCL& matA, VecDescCL& b) const
{
    MLMatrixCL::iterator  itA    = matA.Data.begin();
    MLIdxDescCL::iterator itRow  = matA.RowIdx->begin();
    MLIdxDescCL::iterator itCol  = matA.ColIdx->begin();
    for ( size_t lvl=0; lvl < matA.Data.size(); ++lvl, ++itRow, ++itCol, ++itA)
        SetupSystem_P2( MG_, Coeff_, BndData_, *itA, (lvl == matA.Data.size()-1) ? &b : 0, *itRow, *itCol);
}

//Set up P2-convection
template<class Coeff>
void SetupConvection_P2(const MultiGridCL& MG_, const Coeff& Coeff_, const BndDataCL<> BndData_, MatrixCL& Umat, 
                        VecDescCL* vU, IdxDescCL& RowIdx, IdxDescCL& ColIdx, double tU)
{
  if(vU!=0) vU->Clear(tU);
  MatrixBuilderCL U(&Umat, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
  
  const Uint lvl=RowIdx.TriangLevel();
  const Uint idx=RowIdx.GetIdx();
  
  IdxT Numb[10];
  bool IsOnDirBnd[10];
  
 // Quad3CL<Point3DCL> u;
  LocalP1CL<Point3DCL> Grad[10], GradRef[10];
  SMatrixCL<3,3> T;
  double det;
  double absdet;
  double tmp; //tell the vert points and the edge points
  

  
  for(MultiGridCL::const_TriangTetraIteratorCL sit=MG_.GetTriangTetraBegin(lvl), send=MG_.GetTriangTetraEnd(lvl);
      sit!=send; ++sit)
  {
      GetTrafoTr(T, det, *sit);
      absdet=std::fabs(det);
  
      P2DiscCL::GetGradientsOnRef(GradRef);
      P2DiscCL::GetGradients(Grad, GradRef, T);
       
      Quad3CL<Point3DCL> u(*sit,Coeff_.Vel,tU);
      for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= BndData_.IsOnDirBnd( *sit->GetVertex(i))))
                Numb[i]= sit->GetVertex(i)->Unknowns(idx);
         // u[i]=Coeff_.Vel(sit->GetVertex(i)->GetCoord(),tU);
        }
        //u[4]=Coeff_.Vel(GetBaryCenter(*sit),tU);

        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= BndData_.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }

        //Assemble
        for(int j=0; j<10; ++j)
        {
          const Quad3CL<double> u_Gradj( dot(u, Quad3CL<Point3DCL>(Grad[j], Quad3DataCL::Node))); //??
          if(!IsOnDirBnd[j])
          {
            for(int i=0; i<10; ++i)
            {
              if(!IsOnDirBnd[i])
              {

                 U(Numb[i], Numb[j])+=u_Gradj.quadP2(i, absdet);
              }
            }
          }
          else
              if(vU!=0)
              {
              tmp=j<4? BndData_.GetDirBndValue(*sit->GetVertex(j),tU)
                     : BndData_.GetDirBndValue(*sit->GetEdge(j-4),tU);
                for(int i=0; i<10; ++i)
                {
                if(!IsOnDirBnd[i])
                vU->Data[Numb[i]]-=u_Gradj.quadP2(i,absdet)*tmp;
                }
              }
        }
    
  }
  U.Build();
}


template<class Coeff>
void PoissonP2CL<Coeff>::SetupConvection(MLMatDescCL& matU, VecDescCL& vU, double tU) const
{
  MLMatrixCL ::iterator itU  = matU.Data.begin();
  MLIdxDescCL::iterator itRow= matU.RowIdx->begin();
  MLIdxDescCL::iterator itCol= matU.ColIdx->begin();
  
  for (size_t lvl=0; lvl < matU.Data.size(); ++lvl, ++itRow, ++itCol, ++itU)
     {
     SetupConvection_P2(MG_, Coeff_, BndData_, *itU, (lvl==matU.Data.size()-1) ? &vU:0, *itRow, *itCol, tU);
     }
}

//Initialization of P2 problem
template <class Coeff>
void PoissonP2CL<Coeff>::Init( VecDescCL& vec, instat_scalar_fun_ptr func, double t0) const
{
    Uint lvl= vec.GetLevel(),
         idx= vec.RowIdx->GetIdx();


    for (MultiGridCL::const_TriangVertexIteratorCL sit= const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl), 
        send= const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl); sit != send; ++sit)  //Vertices
    {
        if (sit->Unknowns.Exist(idx))
        {
            vec.Data[sit->Unknowns(idx)]= func( sit->GetCoord(), t0);
        }
    }
    
    for (MultiGridCL::const_TriangEdgeIteratorCL sit= const_cast<const MultiGridCL&>(MG_).GetTriangEdgeBegin(lvl),
        send= const_cast<const MultiGridCL&>(MG_).GetTriangEdgeEnd(lvl); sit!=send; ++sit)      //Edges
    {
        if ( sit->Unknowns.Exist(idx))
            vec.Data[sit->Unknowns(idx)]= func( GetBaryCenter( *sit), t0);
    }

}

//Setup Rhs of instat P2
template<class Coeff>
void PoissonP2CL<Coeff>::SetupInstatRhs(VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf) const
{
  vA.Clear(tA);
  vM.Clear(tA);
  vf.Clear(tf);
  
  const Uint lvl = vA.GetLevel(),
             idx = vA.RowIdx->GetIdx();
             
  Comment("InstatPoissonP2CL::SetupInstatRhs with Index"<<idx<<std::endl, DebugNumericC);
  Quad2CL<> rhs;
  IdxT Numb[10];
  bool IsOnDirBnd[10];


  SMatrixCL<3, 5> Grad[10], GradRef[10];
  SMatrixCL<3, 3> T;
  double coup[10][10];
  double det, absdet;
        
    //fill value part of matrices
    GetGradientsOnRef(GradRef);
    for(MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), 
                            send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl); sit != send; ++sit)
    {
     
        GetTrafoTr(T, det, *sit);
        MakeGradients(Grad, GradRef, T);
        absdet = std:: fabs(det);
        
        //collet some information about the edges and verts of the tetra and save it in the Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= BndData_.IsOnDirBnd( *sit->GetVertex(i))))
                Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        }

        for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= BndData_.IsOnDirBnd( *sit->GetEdge(i))))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }
        //compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                coup[i][j] = Coeff::alpha * QuadGrad(Grad, i, j)*absdet;
                coup[i][j]+= Quad(*sit, Coeff::q, i, j, 0.0)*absdet;
                coup[j][i] = coup[i][j];
            }
        //assemble
        for(int j=0; j<10; ++j)
            if(IsOnDirBnd[j])    //Vert/Edge j is on a Dirichlet Boundary
            {//coupling with vertex j on right-hand-side
              double tmp;
              tmp=j<4? BndData_.GetDirBndValue(*sit->GetVertex(j),tA)
                     : BndData_.GetDirBndValue(*sit->GetEdge(j-4),tA);
                for(int i=0; i<10; ++i)
                {
                    if(!IsOnDirBnd[i]) //Vert/Edge j is not on a Dirichlet Boundary
                    {
                       vA.Data[Numb[i]]-=coup[j][i]*tmp;
                       vM.Data[Numb[i]]-=P2DiscCL::GetMass(i, j)*absdet*tmp;
                    }
                }
            }
        //rhs.assign(*sit, Coeff_.f, tf);

        for(int i=0; i<10;++i)
        {    
            if(!IsOnDirBnd[i])
            {
                double valf;
                valf=Quad(*sit, Coeff::f, i, tf)*absdet;
                vf.Data[Numb[i]]+=valf;                       //rhs.quadP2(i, absdet);
            }
            if ( i<4 ? BndData_.IsOnNatBnd(*sit->GetVertex(i))
                   : BndData_.IsOnNatBnd(*sit->GetEdge(i-4)) )
            {
              Uint face;
              int n=i<4? 3:2;  //Three faces for a vertex and two faces for a point on the edge;
              for (int f=0; f < n; ++f)
                {
                 face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);  
                 if ( sit->IsBndSeg(face))
                    {
                     vA.Data[Numb[i]]+=QuadP2(*sit, face, BndData_, i);
                    }
                }
            }
            
        }
    }
}
//Setup Instat P2 problem
template<class Coeff>
void SetupInstatSystem_P2( const MultiGridCL& MG, const Coeff&, const BndDataCL<> BndData_, MatrixCL& Amat, MatrixCL& Mmat, IdxDescCL& RowIdx, IdxDescCL& ColIdx)
{ //Setup stiffness and mass matrices
    
    const Uint lvl = RowIdx.TriangLevel(),
               idx = RowIdx.GetIdx();
    MatrixBuilderCL A(&Amat, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    MatrixBuilderCL M(&Mmat, RowIdx.NumUnknowns(), ColIdx.NumUnknowns());
    
    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
#ifndef _PAR
 std::cout<< "entering SetupInstatSystem: " << ColIdx.NumUnknowns()<< " unknowns,"<<std::endl;
#endif

//fill value part of matrices
    SMatrixCL<3, 5> Grad[10], GradRef[10];
    SMatrixCL<3, 3> T;
    double coup[10][10];
    double det, absdet;
    
    GetGradientsOnRef(GradRef);
    for(MultiGridCL::const_TriangTetraIteratorCL sit=MG.GetTriangTetraBegin(lvl), send=MG.GetTriangTetraEnd(lvl);
    sit != send; ++sit)
    {
        
        GetTrafoTr(T,det, *sit);
        MakeGradients(Grad, GradRef, T);
        absdet = std:: fabs(det);
        //collet some information about the edges and verts of the tetra and save it in the Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if (!(IsOnDirBnd[i]=BndData_.IsOnDirBnd( *sit-> GetVertex(i))))
                Numb[i]=sit->GetVertex(i)->Unknowns(idx);
        }
        for(int i=0; i<6; ++i)
        {
            if(!(IsOnDirBnd[i+4]=BndData_.IsOnDirBnd( *sit-> GetEdge(i))))
                Numb[i+4]=sit->GetEdge(i)->Unknowns(idx);
        }
        //compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                coup[i][j] = Coeff::alpha * QuadGrad(Grad, i, j)*absdet;
                coup[i][j]+= Quad(*sit, Coeff::q, i, j, 0.0)*absdet;  //reaction term
                coup[j][i] = coup[i][j];
            }
        //assemble
        for(int i=0; i<10; ++i)
            if(!IsOnDirBnd[i])    //Vert/Edge i is not on a Dirichlet Boundary
                for(int j=0; j<10; ++j)
                {
                    if(!IsOnDirBnd[j]) //Vert/Edge j is not on a Dirichlet Boundary
                    {
                        A(Numb[i], Numb[j])+=coup[j][i];
                        M(Numb[i], Numb[j])+=P2DiscCL::GetMass(i, j)*absdet;
                    }
                }
                
    }
#ifndef _PAR
    std::cout<< "done: value part fill"<<std::endl;
#endif

   A.Build();
   M.Build();

}

template<class Coeff>
void PoissonP2CL<Coeff>::SetupInstatSystem( MLMatDescCL& matA, MLMatDescCL& matM, double) const
{
    MLIdxDescCL::iterator itRow  = matA.RowIdx->begin();
    MLIdxDescCL::iterator itCol  = matA.ColIdx->begin();
    MLMatrixCL::iterator  itM    = matM.Data.begin();
    for ( MLMatrixCL::iterator itA= matA.Data.begin(); itA != matA.Data.end(); ++itA, ++itM, ++itRow, ++itCol)
        SetupInstatSystem_P2( MG_, Coeff_, BndData_, *itA, *itM, *itRow, *itCol);
}

/// \todo CheckSolution checks 2-norm and max-norm just on vertices and not on edges
template <class Coeff>
double PoissonP2CL<Coeff>::CheckSolution(const VecDescCL& lsg, instat_scalar_fun_ptr Lsg, double t) const
{
  double diff, maxdiff=0, norm2= 0, L2=0;
  Uint lvl=lsg.GetLevel(),
       Idx=lsg.RowIdx->GetIdx();

  const_DiscSolCL sol(&lsg, &GetBndData(), &GetMG());

  std::cout << "Difference to exact solution:" << std::endl;

  for (MultiGridCL::const_TriangTetraIteratorCL
    sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl),
    send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
    sit != send; ++sit)
  {
    double absdet= sit->GetVolume()*6.,
           sum= 0;

    for(Uint i=0; i<4; ++i)
    {
      diff= (sol.val(*sit->GetVertex(i)) - Lsg(sit->GetVertex(i)->GetCoord(),t));
      sum+= diff*diff;
    }
    sum/= 120;
    diff= sol.val(*sit, 0.25, 0.25, 0.25) - Lsg(GetBaryCenter(*sit),t);
    sum+= 2./15. * diff*diff;
    L2+= sum*absdet;
  }
#ifdef _PAR
  L2= ProcCL::GlobalSum(L2);
#endif

  L2= std::sqrt(L2);

  for (MultiGridCL::const_TriangVertexIteratorCL
    sit=const_cast<const MultiGridCL&>(MG_).GetTriangVertexBegin(lvl),
    send=const_cast<const MultiGridCL&>(MG_).GetTriangVertexEnd(lvl);
    sit != send; ++sit)
  {
    if (sit->Unknowns.Exist(Idx))
#ifdef _PAR
      if (sit->IsExclusive())
#endif
      {
        diff= std::fabs( Lsg(sit->GetCoord(),t) - lsg.Data[sit->Unknowns(Idx)] );
        norm2+= diff*diff;
        if (diff>maxdiff)
        {
          maxdiff= diff;
        }
      }
  }

  int Lsize = lsg.Data.size();
#ifdef _PAR
  Lsize   = lsg.RowIdx->GetGlobalNumUnknowns(MG_);
  norm2   = ProcCL::GlobalSum(norm2, Drops_MasterC);
  maxdiff = ProcCL::GlobalMax(maxdiff, Drops_MasterC);
#endif
  std::cout << "  2-Norm= " << std::sqrt(norm2)
            << "\nw-2-Norm= " << std::sqrt(norm2/Lsize)
            << "\nmax-Norm= " << maxdiff
            << "\n L2-Norm= " << L2 << std::endl;
  return L2;
}

template<class Coeff>
void PoissonP2CL<Coeff>::SetNumLvl( size_t n)
{
#ifdef _PAR
    if (n>1)
        throw DROPSErrCL(" PoissonP2CL<Coeff>::SetNumLvl: Sorry, in parallel version no multi-level is implemented, yet!");
#endif
    match_fun match= MG_.GetBnd().GetMatchFun();
    idx.resize( n, P2_FE, BndData_, match);
    A.Data.resize( idx.size());
}

//========================================================================================================
//
//                           Marker classes for adaptive refinement
//
//========================================================================================================


template <class _TetraEst, class _ProblemCL>
template <class BndData_, class _VD>
void PoissonErrEstCL<_TetraEst, _ProblemCL>::Init(const P1EvalCL<double, BndData_, _VD>& sol)
{
    MultiGridCL& mg= _Problem.GetMG();
    const Uint lvl= sol.GetLevel();
    const typename _ProblemCL::BndDataCL& bnd= _Problem.GetBndData();

    double tmp;
    _InitGlobErr= 0;
    for (MultiGridCL::TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        tmp= _Estimator( *sit, *sol.GetSolution(), bnd);
        _InitGlobErr+= tmp*tmp;
    }
    _InitGlobErr= std::sqrt(_InitGlobErr);
    _ActGlobErr= _InitGlobErr;
    if (_outp)
        *_outp << "ErrorEstimator is initialized now; initial error: " << _InitGlobErr
               << " We try to reduce by a factor of " << _RelReduction
               << "." << std::endl;
}

template <class _TetraEst, class _ProblemCL>
template <class BndData_, class _VD>
bool PoissonErrEstCL<_TetraEst, _ProblemCL>::Estimate(const P1EvalCL<double, BndData_, const _VD>& sol)
{
    const MultiGridCL& mg= _Problem.GetMG();
    const typename _ProblemCL::BndDataCL& bnd= _Problem.GetBndData();
    const Uint lvl= sol.GetLevel();
    const VecDescCL& lsg= *sol.GetSolution();
    const double exp_err= _InitGlobErr*_RelReduction/std::sqrt(_Meas);
    const double unref_bnd= exp_err/2./std::pow(2, _ConvExponent);
    double globalerr= 0;
    double localerr;
    double localerr_dist;
    Uint num_tetra= 0;

    _NumLastMarkedForRef= 0;
    _NumLastMarkedForDel= 0;
    for (MultiGridCL::const_TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        ++num_tetra;
        localerr= _Estimator(*sit, lsg, bnd);
        localerr_dist= localerr/std::sqrt(sit->GetVolume());
        globalerr+= localerr*localerr;
        if ( sit->IsUnrefined() ) {
            if (localerr_dist>exp_err)
            {
                sit->SetRegRefMark();
                ++_NumLastMarkedForRef;
            }
            else if (localerr_dist<unref_bnd)
                {
                    sit->SetRemoveMark();
                    ++_NumLastMarkedForDel;
                }
        }
    }
    globalerr= std::sqrt(globalerr);
    if (globalerr<_InitGlobErr*_RelReduction)
    {
        for (MultiGridCL::const_TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
             sit != send; ++sit)
            if (sit->IsUnrefined()) sit->SetNoRefMark();
        _NumLastMarkedForDel= _NumLastMarkedForRef= 0;
    }
    if (_outp)
        *_outp << "Estimated global W1_2-error: " << globalerr << ". This is a reduction of " << globalerr/_ActGlobErr
               << " compared to the last estimation. Marked " << _NumLastMarkedForRef
               << ", unmarked " << _NumLastMarkedForDel << " out of " << num_tetra << " tetrahedrons."
               << std::endl;
    _ActGlobErr= globalerr;
    return _NumLastMarkedForRef || _NumLastMarkedForDel;
}

template <class _TetraEst, class _ProblemCL>
template <class BndData_, class _VD>
void DoerflerMarkCL<_TetraEst, _ProblemCL>::Init(const P1EvalCL<double, BndData_, _VD>& sol)
{
    MultiGridCL& mg= _Problem.GetMG();
    const Uint lvl= sol.GetLevel();
    const typename _ProblemCL::BndDataCL& bnd= _Problem.GetBndData();

    _InitGlobErr= 0;
    for (MultiGridCL::TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        _InitGlobErr+= _Estimator( *sit, *sol.GetSolution(), bnd);
    }
    _InitGlobErr= std::sqrt(_InitGlobErr);
    _ActGlobErr= _InitGlobErr;
    if (_outp)
        *_outp << "ErrorEstimator is initialized now; initial error: " << _InitGlobErr
               << " We try to reduce by a factor of " << _RelReduction
               << " by marking the tetrahedrons with largest error-estimates, until the error marked"
               << " is at least " << _Threshold << " of the actual global error." << std::endl;
}

template <class _TetraEst, class _ProblemCL>
template <class BndData_, class _VD>
bool DoerflerMarkCL<_TetraEst, _ProblemCL>::Estimate(const P1EvalCL<double, BndData_, const _VD>& sol)
{
    Err_ContCL err_est;
    const MultiGridCL& mg= _Problem.GetMG();
    const Uint lvl= sol.GetLevel();
    const VecDescCL& lsg= *sol.GetSolution();
    const typename _ProblemCL::BndDataCL& bnd= _Problem.GetBndData();

    _NumLastMarkedForRef= 0;
    _NumLastMarkedForDel= 0;
    for (MultiGridCL::const_TriangTetraIteratorCL sit=mg.GetTriangTetraBegin(lvl), send=mg.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        const double localerr= _Estimator(*sit, lsg, bnd);
        err_est.push_back( std::make_pair(&*sit, localerr) );
    }
    const double globalerr_sq= std::accumulate(err_est.begin(), err_est.end(), 0.0, AccErrCL() );
    const double globalerr= std::sqrt(globalerr_sq);
    const double ref_threshold_sq= globalerr_sq*_Threshold*_Threshold;
    if (globalerr>=_InitGlobErr*_RelReduction && _DoMark)
    {
        std::sort( err_est.begin(), err_est.end(), Err_Pair_GTCL() );
        double akt_ref_err_sq= 0;
        const Uint min_tetra= static_cast<Uint>(err_est.size()*_min_tetra_ratio);
        for (Err_ContCL::iterator it= err_est.begin(), theend= err_est.end();
             it != theend && (akt_ref_err_sq < ref_threshold_sq || _NumLastMarkedForRef < min_tetra); ++it)
            if ( it->first->IsUnrefined() )
            {
                it->first->SetRegRefMark();
                akt_ref_err_sq+= it->second;
                ++_NumLastMarkedForRef;
            }
    }
    if (_outp)
        *_outp << "Estimated global W1_2-error: " << globalerr << ". This is a reduction of " << globalerr/_ActGlobErr
               << " compared to the last estimation. Marked " << _NumLastMarkedForRef
               << ", unmarked " << _NumLastMarkedForDel << " out of " << err_est.size() << " tetrahedrons, "
               << "which account for " << _Threshold << " of the global error."
               << std::endl;
    _ActGlobErr= globalerr;
    return _NumLastMarkedForRef || _NumLastMarkedForDel;
}


} // end of namespace DROPS

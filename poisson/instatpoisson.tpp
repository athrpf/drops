//**************************************************************************
// File:    instatpoisson.tpp                                              *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - November, 11 2002                                      *
//**************************************************************************

#include "num/discretize.h"

namespace DROPS
{


template<class MGB, class Coeff>
void InstatPoissonP1CL<MGB,Coeff>::CreateNumbering(Uint level, IdxDescCL* idx)
// used for numbering of the Unknowns depending on the index IdxDesc[idxnum].
// sets up the description of the index idxnum in IdxDesc[idxnum],
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, that IdxDesc[idxnr].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx->TriangLevel = level;
    idx->NumUnknowns = 0;

    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL

    // allocate space for indices; number unknowns in TriangLevel level
    CreateNumbOnVertex( idxnum, idx->NumUnknowns, idx->NumUnknownsVertex,
                        _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level), GetBndData() );
}


template<class MGB, class Coeff>
void InstatPoissonP1CL<MGB,Coeff>::DeleteNumbering(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level) );
}

//========================================================
//
//                Set up matrices and rhs
//
//========================================================

inline double Quad2D(const TetraCL& t, Uint face, Uint vert, InstatPoissonBndDataCL::bnd_val_fun bfun, double time)
// Integrate neu_val() * phi_vert over face
{
    const BndIdxT bidx= t.GetBndIdx(face);
    Point2DCL vc2D[3];
    const VertexCL* v[3];
    const BndPointSegEqCL comp(bidx);
    
    v[0]= t.GetVertex(vert);
    for (Uint i=0, k=1; i<3; ++i)
    {
        if (VertOfFace(face,i)!=vert)
            v[k++]= t.GetVertex( VertOfFace(face,i) );
        vc2D[i]= std::find_if( v[i]->GetBndVertBegin(), v[i]->GetBndVertEnd(), comp)->GetCoord2D();
    }
    const double f0= bfun(vc2D[0], time);
    const double f1= bfun(vc2D[1], time) +  bfun( vc2D[2], time);
    const double f2= bfun(1./3.*(vc2D[0] + vc2D[1] + vc2D[2]), time);    //Barycenter of Face
    const double absdet= FuncDet2D(v[1]->GetCoord() - v[0]->GetCoord(), v[2]->GetCoord() - v[0]->GetCoord());
    return (11./240.*f0 + 1./240.*f1 + 9./80.*f2) * absdet;
}


template<class MGB, class Coeff>
void InstatPoissonP1CL<MGB,Coeff>::SetupInstatRhs(VecDescCL& vA, VecDescCL& vM, double tA, VecDescCL& vf, double tf) const
// Sets up the time dependent right hand sides including couplings 
// resulting from inhomogeneous dirichlet bnd conditions
{
  vA.Clear();
  vM.Clear();
  vf.Clear();
    
  const Uint lvl = vA.RowIdx->TriangLevel,
             idx = vA.RowIdx->GetIdx();
  SMatrixCL<3,4> G;

  double coup[4][4];
  double det;
  double absdet;
  IdxT UnknownIdx[4];
	
	StripTimeCL strip( &_Coeff.f, tf);

  for (MultiGridCL::const_TriangTetraIteratorCL
    sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl),
    send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
    sit != send; ++sit)
  {
    P1DiscCL::GetGradients(G,det,*sit);
    absdet= fabs(det);

    for(int i=0; i<4; ++i)
    {
      for(int j=0; j<=i; ++j)
      {
        // dot-product of the gradients
        coup[i][j]= ( G(0,i)*G(0,j)+G(1,i)*G(1,j)+G(2,i)*G(2,j) )/6.0*absdet;
        // coup[i][j]+= P1DiscCL::Quad(*sit, &_Coeff.q, i, j)*absdet;
        coup[j][i]= coup[i][j];
      }
      UnknownIdx[i]= sit->GetVertex(i)->Unknowns.Exist() ? sit->GetVertex(i)->Unknowns(idx) 
                                                             : -1ul;
    }
    for(int i=0; i<4; ++i)    // assemble row i
      if (sit->GetVertex(i)->Unknowns.Exist())  // vertex i is not on a Dirichlet boundary
      {
        for(int j=0; j<4;++j)
        {
          if (!sit->GetVertex(j)->Unknowns.Exist()) // vertex j is on a Dirichlet boundary
          { 
						// coupling with vertex j on right-hand-side
            vA.Data[UnknownIdx[i]]-= coup[j][i] * _BndData.GetDirBndValue(*sit->GetVertex(j), tA);
            vM.Data[UnknownIdx[i]]-= (i==j ? 1./60. : 1./120.)*absdet 
                                     * _BndData.GetDirBndValue(*sit->GetVertex(j), tA);
          }
        }
        
        vf.Data[UnknownIdx[i]]+= P1DiscCL::Quad(*sit, &strip.GetFunc, i)*absdet;
        if ( _BndData.IsOnNeuBnd(*sit->GetVertex(i)) )
          for (int f=0; f < 3; ++f)
            if ( sit->IsBndSeg(FaceOfVert(i, f)) )
              vA.Data[UnknownIdx[i]]+=
                Quad2D(*sit, FaceOfVert(i, f), i, _BndData.GetBndFun( sit->GetBndIdx( FaceOfVert(i,f))), tA);
      }
  }
}


template<class MGB, class Coeff>
void InstatPoissonP1CL<MGB,Coeff>::SetupInstatSystem( MatDescCL& Amat, MatDescCL& Mmat) const
// Sets up the stiffness matrix
{
  MatrixBuilderCL A( &Amat.Data, Amat.RowIdx->NumUnknowns, Amat.ColIdx->NumUnknowns);
	MatrixBuilderCL M( &Mmat.Data, Mmat.RowIdx->NumUnknowns, Mmat.ColIdx->NumUnknowns);
	 
  const Uint lvl = Amat.RowIdx->TriangLevel;
  const Uint idx = Amat.RowIdx->GetIdx();

  SMatrixCL<3,4> G;
    
  double coup[4][4];
  double det;
  double absdet;
  IdxT UnknownIdx[4];
	 
  for (MultiGridCL::const_TriangTetraIteratorCL
    sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), 
    send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
    sit != send; ++sit)
  {
    P1DiscCL::GetGradients(G,det,*sit);
    absdet= fabs(det);

    for(int i=0; i<4; ++i)
    {
      for(int j=0; j<=i; ++j)
      {
        // dot-product of the gradients
        coup[i][j]= ( G(0,i)*G(0,j)+G(1,i)*G(1,j)+G(2,i)*G(2,j) )/6.0*absdet;
        // coup[i][j]+= P1DiscCL::Quad(*sit, &_Coeff.q, i, j)*absdet;
        coup[j][i]= coup[i][j];
      }
      UnknownIdx[i]= sit->GetVertex(i)->Unknowns.Exist() ? sit->GetVertex(i)->Unknowns(idx) : -1ul;
    }
		  
    for(int i=0; i<4; ++i)    // assemble row i
      if (sit->GetVertex(i)->Unknowns.Exist())  // vertex i is not on a Dirichlet boundary
      {
        for(int j=0; j<4;++j)
        {
          if (sit->GetVertex(j)->Unknowns.Exist()) // vertex j is not on a Dirichlet boundary
          {
            A( UnknownIdx[i], UnknownIdx[j])+= coup[j][i]; 
            M( UnknownIdx[i], UnknownIdx[j])+= (i==j ? 1./60. : 1./120.)*absdet;
          }
          // else coupling with vertex j on right-hand-side  --> 0
        }
      }
  }
  A.Build();
  M.Build();
}


template<class MGB, class Coeff>
void InstatPoissonP1CL<MGB,Coeff>::Init( VecDescCL& vec, scalar_instat_fun_ptr func, double t0) const
{
    Uint lvl= vec.RowIdx->TriangLevel,
         idx= vec.RowIdx->GetIdx();
    
    
    for (MultiGridCL::const_TriangVertexIteratorCL sit= _MG.GetTriangVertexBegin(lvl), send= _MG.GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (sit->Unknowns.Exist(idx))
        {
            vec.Data[sit->Unknowns(idx)]= func( sit->GetCoord(), t0);
        }
    }
    
}

//========================================================
//
//                Check solution 
//
//========================================================

template<class MGB, class Coeff>
void InstatPoissonP1CL<MGB,Coeff>::CheckSolution(const VecDescCL& lsg, 
  scalar_instat_fun_ptr Lsg, double t) const
{
  double diff, maxdiff=0, norm2= 0, L2=0;
  Uint lvl=lsg.RowIdx->TriangLevel,
       Idx=lsg.RowIdx->GetIdx();
  
  DiscSolCL sol(&lsg, &GetBndData(), &GetMG(), t);
    
  std::cerr << "Abweichung von der tatsaechlichen Loesung:" << std::endl;

  for (MultiGridCL::const_TriangTetraIteratorCL 
    sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl),
    send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
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
  L2= sqrt(L2);

  for (MultiGridCL::const_TriangVertexIteratorCL 
    sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl),
    send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
    sit != send; ++sit)
  {
    if (sit->Unknowns.Exist())
    {
      diff= fabs( Lsg(sit->GetCoord(),t) - lsg.Data[sit->Unknowns(Idx)] );
      norm2+= diff*diff;
      if (diff>maxdiff)
      {
        maxdiff= diff;
      }
    }
  }
  std::cerr << "  2-Norm= " << ::sqrt(norm2)                 << std::endl
            << "w-2-Norm= " << ::sqrt(norm2/lsg.Data.size()) << std::endl
            << "max-Norm= " << maxdiff                       << std::endl
            << " L2-Norm= " << L2                            << std::endl;

}


} // end of namespace DROPS

//**************************************************************************
// File:    poisson.tpp                                                    *
// Content: classes that constitute the poisson-problem                    *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
// Version: 0.1                                                            *
// History: begin - March, 16 2001                                         *
//**************************************************************************

#include "num/discretize.h"

namespace DROPS
{

template<class Coeff>
void PoissonP1CL<Coeff>::CreateNumbering(Uint level, IdxDescCL* idx)
// used for numbering of the Unknowns depending on the index IdxDesc[idxnum].
// sets up the description of the index idxnum in IdxDesc[idxnum],
// allocates memory for the Unknown-Indices on TriangLevel level und numbers them.
// Remark: expects, thatIdxDesc[idxnr].NumUnknownsVertex etc. are set.
{
    // set up the index description
    idx->TriangLevel = level;
    idx->NumUnknowns = 0;

    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL

    // allocate space for indices; number unknowns in TriangLevel level
    CreateNumbOnVertex( idxnum, idx->NumUnknowns, idx->NumUnknownsVertex,
                        _MG.GetTriangVertexBegin( level), _MG.GetTriangVertexEnd( level), GetBndData());
}


template<class Coeff>
void PoissonP1CL<Coeff>::DeleteNumbering(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetTriangVertexBegin( level), _MG.GetTriangVertexEnd( level) );
}

//========================================================
//
//                Set up matrices and rhs
//
//========================================================


inline double Quad2D(const TetraCL& t, Uint face, Uint vert, PoissonBndDataCL::bnd_val_fun bfun)
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
        //std::cerr << "\nVertex " << i << ": " << vc2D[i] << std::endl;
    }
    const double f0= bfun(vc2D[0]);
    const double f1= bfun(vc2D[1]) +  bfun( vc2D[2]);
    const double f2= bfun(1./3.*(vc2D[0] + vc2D[1] + vc2D[2]) );    //Barycenter of Face
    const double absdet= FuncDet2D(v[1]->GetCoord() - v[0]->GetCoord(), v[2]->GetCoord() - v[0]->GetCoord());
    //std::cerr << '\n' << f0 << ' ' << f1 << ' ' << f2 << ' ' << absdet << std::endl;
    //std::cerr << '\n' << v[0]->GetCoord() << ' ' << v[1]->GetCoord() << ' ' << v[2]->GetCoord() << ' ' << absdet << std::endl;
    return (11./240.*f0 + 1./240.*f1 + 9./80.*f2) * absdet;
}

/*
inline double Quad2D(const TetraCL& t, Uint face, Uint vert, PoissonBndDataCL::bnd_val_fun bfun)
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
        //std::cerr << "\nVertex3D " << i << ": " << v[i]->GetCoord();
        //std::cerr << "\nVertex2D " << i << ": " << vc2D[i];
    }
    const double f0= bfun(vc2D[0]);
    const double f1= bfun(vc2D[1]);
    const double f2= bfun(vc2D[2]);    
    const double absdet= FuncDet2D(v[1]->GetCoord() - v[0]->GetCoord(), v[2]->GetCoord() - v[0]->GetCoord());
    //std::cerr << '\n' << f0 << ' ' << f1 << ' ' << f2 << ' ' << std::endl;
    //std::cerr << '\n' << v[0]->GetCoord() << ' ' << v[1]->GetCoord() << ' ' << v[2]->GetCoord() << ' ' << absdet << std::endl;
    //std::cerr << absdet <<"\n";
    //std::cerr << (1./12.*f0 + 1./24.*f1 + 1./24.*f2) << "\n";
    //std::cerr << (1./12.*f0 + 1./24.*f1 + 1./24.*f2)*absdet << "\n";
    return (1./12.*f0 + 1./24.*f1 + 1./24.*f2) * absdet;
}
*/

template<class Coeff>
void PoissonP1CL<Coeff>::SetupSystem(MatDescCL& Amat, VecDescCL& b) const
// Sets up the stiffness matrix and right hand side
{
    b.Clear();
    
    const Uint lvl    = Amat.GetRowLevel(),
               idx    = Amat.RowIdx->GetIdx();
    MatrixBuilderCL A( &Amat.Data, Amat.RowIdx->NumUnknowns, Amat.ColIdx->NumUnknowns); 

//    SMatrixCL<3,3> T;
//    SMatrixCL<3,4> Gref(0.0); 
//    Gref(0,1)= Gref(1,2)= Gref(2,3)= 1.; // gradients on ref. tetra
//    Gref(0,0)= Gref(1,0)= Gref(2,0)= -1.;
    SMatrixCL<3,4> G;

    double coup[4][4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        P1DiscCL::GetGradients(G,det,*sit);   
//        SMatrixCL<3,4> G= T*Gref;
/*        G(0,0)= -T(0,0)-T(0,1)-T(0,2);
        G(0,1)= T(0,0); G(0,2)= T(0,1); G(0,3)= T(0,2);
        G(1,0)= -T(1,0)-T(1,1)-T(1,2);
        G(1,1)= T(1,0); G(1,2)= T(1,1); G(1,3)= T(1,2);
        G(2,0)= -T(2,0)-T(2,1)-T(2,2);
        G(2,1)= T(2,0); G(2,2)= T(2,1); G(2,3)= T(2,2);
*/        absdet= fabs(det);

        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= ( G(0,i)*G(0,j) + G(1,i)*G(1,j) + G(2,i)*G(2,j) )/6.0*absdet;
                coup[i][j]+= P1DiscCL::Quad(*sit, &_Coeff.q, i, j)*absdet;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetVertex(i)->Unknowns.Exist(idx) ? sit->GetVertex(i)->Unknowns(idx) 
                                                                  : NoIdx;
        }
        for(int i=0; i<4; ++i)    // assemble row i
            if (sit->GetVertex(i)->Unknowns.Exist(idx))  // vertex i is not on a Dirichlet boundary
            {
                for(int j=0; j<4;++j)
                {
                    if (sit->GetVertex(j)->Unknowns.Exist(idx)) // vertex j is not on a Dirichlet boundary
                    {
                        A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i]; 
                    }
                    else // coupling with vertex j on right-hand-side
                    {
                        b.Data[UnknownIdx[i]]-= coup[j][i] * _BndData.GetDirBndValue(*sit->GetVertex(j));
                    }
                }
                b.Data[UnknownIdx[i]]+= P1DiscCL::Quad(*sit, &_Coeff.f, i)*absdet;
                if ( _BndData.IsOnNeuBnd(*sit->GetVertex(i)) )
                    for (int f=0; f < 3; ++f)
                        if ( sit->IsBndSeg(FaceOfVert(i, f)) )
                            b.Data[UnknownIdx[i]]+= Quad2D(*sit, FaceOfVert(i, f), i, _BndData.GetBndFun(sit->GetBndIdx(FaceOfVert(i,f))) );
            }
    }
    
    A.Build();
}



template<class Coeff>
void PoissonP1CL<Coeff>::SetupStiffnessMatrix(MatDescCL& Amat) const
// Sets up the stiffness matrix
{
    MatrixBuilderCL A(&Amat.Data, Amat.RowIdx->NumUnknowns, Amat.ColIdx->NumUnknowns);
    const Uint lvl    = Amat.GetRowLevel();
    const Uint idx    = Amat.RowIdx->GetIdx();

    SMatrixCL<3,4> G;
    
    double coup[4][4];
    double det;
    double absdet;
    IdxT UnknownIdx[4];

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        P1DiscCL::GetGradients(G,det,*sit);
        absdet= fabs(det);
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<=i; ++j)
            {
                // dot-product of the gradients
                coup[i][j]= ( G(0,i)*G(0,j) + G(1,i)*G(1,j) + G(2,i)*G(2,j) )/6.0*absdet;
                coup[i][j]+= P1DiscCL::Quad(*sit, &_Coeff.q, i, j)*absdet;
                coup[j][i]= coup[i][j];
            }
            UnknownIdx[i]= sit->GetVertex(i)->Unknowns.Exist(idx) ? sit->GetVertex(i)->Unknowns(idx) : NoIdx;
        }
        for(int i=0; i<4; ++i)    // assemble row i
            if (sit->GetVertex(i)->Unknowns.Exist(idx))  // vertex i is not on a Dirichlet boundary
            {
                for(int j=0; j<4;++j)
                {
                    if (sit->GetVertex(j)->Unknowns.Exist(idx)) // vertex j is not on a Dirichlet boundary
                    {
                        A(UnknownIdx[i], UnknownIdx[j])+= coup[j][i]; 
                    }
                    // else coupling with vertex j on right-hand-side  --> 0
                }
            }
    }
    A.Build();
}

template<class Coeff>
void PoissonP1CL<Coeff>::SetupProlongation(MatDescCL& P, IdxDescCL* cIdx, IdxDescCL* fIdx) const
// This only works, if Interpolate is called after every refinement of the multigrid.
{
    SetupP1ProlongationMatrix( _MG, P, cIdx, fIdx);
}

//===================================================
//
//   check computed solution, estimate error etc.
//
//===================================================

template<class Coeff>
double PoissonP1CL<Coeff>::CheckSolution(const VecDescCL& lsg, scalar_fun_ptr Lsg) const
{
    double diff, maxdiff=0, norm2= 0, L2=0;
    Uint lvl=lsg.GetLevel(),
         Idx=lsg.RowIdx->GetIdx();
    
    DiscSolCL sol(&lsg, &GetBndData(), &GetMG());
    
    std::cerr << "Abweichung von der tatsaechlichen Loesung:" << std::endl;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double absdet= sit->GetVolume()*6.,
            sum= 0;
        for(Uint i=0; i<4; ++i)
        {
            diff= (sol.val(*sit->GetVertex(i)) - Lsg(sit->GetVertex(i)->GetCoord()));
            sum+= diff*diff;
        }
        sum/= 120;
        diff= sol.val(*sit, 0.25, 0.25, 0.25) - Lsg(GetBaryCenter(*sit));
        sum+= 2./15. * diff*diff;
        L2+= sum*absdet;
    }
    L2= sqrt(L2);

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (sit->Unknowns.Exist(Idx))
        {
           diff= fabs( Lsg(sit->GetCoord()) - lsg.Data[sit->Unknowns(Idx)] );
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
    return L2;
}

template<class Coeff>
void PoissonP1CL<Coeff>::GetDiscError(const MatDescCL& A, scalar_fun_ptr Lsg) const
{
    Uint lvl= A.GetColLevel(),
         idx= A.ColIdx->GetIdx();
    VectorCL lsg(A.Data.num_cols());

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (sit->Unknowns.Exist(idx))
        {
            lsg[sit->Unknowns(idx)]= Lsg(sit->GetCoord());
        }
    }
    
    VectorCL res= A.Data*lsg-b.Data; 
    std::cerr <<"|| Ax - b || = "<< res.norm()<<", max "<<res.supnorm()
              <<" (with x=continuous solution)"<<std::endl;
}



template<class Coeff>
bool PoissonP1CL<Coeff>::EstimateError (const VecDescCL& lsg, const double rel_loc_tol, double& globalerr, est_fun est)
{
    const Uint lvl= lsg.GetLevel();
    Uint num_ref= 0;
    Uint num_un_ref= 0;
    const Uint num_tetra= std::distance(_MG.GetTriangTetraBegin(lvl), _MG.GetTriangTetraEnd(lvl));
    const double exp_err= globalerr/0.875; // divident is the volume of the domain
    double localerr;
    double localerr_dist;

    globalerr= 0.0;
    for (MultiGridCL::TriangTetraIteratorCL sit=_MG.GetTriangTetraBegin(lvl), send=_MG.GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        localerr= est(*sit, lsg, GetBndData());
        localerr_dist= localerr/::sqrt(sit->GetVolume());
        globalerr+= localerr*localerr;
        if ( sit->IsUnrefined() )
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
    globalerr= ::sqrt(globalerr);
    std::cerr << "Estimated global W1_2-error: " << globalerr << ". Marked " << num_ref << ", unmarked " << num_un_ref
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
    const double absdet= fabs(det);
    // Integral over tetra
    circumcircle(t, cc_center, cc_radius);
    const double P0f= Quad3CL::Quad(t, &CoeffCL::f)*6.; // absdet/vol = 6.
    _err+= 4.*cc_radius*cc_radius*P0f*P0f*absdet/6.;
//    _err+= 4.*cc_radius*cc_radius*P1DiscCL::norm_L2_sq(t, &CoeffCL::f)*absdet;

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
            const typename BndDataCL::BndSegDataCL* bdat= &Bnd.GetSegData(bidx);
            if ( bdat->IsNeumann() )
            {
//std::cerr << "neumann";
                Point2DCL vc2D[3];
                const VertexCL* v[4];

                const BndPointSegEqCL comp(bidx);
                double n_du= 0.0;  // scalar product of Du and the normal
                for (int i=0; i<4; ++i)
                {
                    v[i]= t.GetVertex(i);
                    n_du+= (v[i]->Unknowns.Exist(Idx) ? sol.Data[v[i]->Unknowns(Idx)]
                               : Bnd.GetDirBndValue(*v[i]) )
                          *(H(0,i)*normal[0] + H(1,i)*normal[1] + H(2,i)*normal[2]);
                }
                for (int i=0; i<3; ++i)
                    vc2D[i]= std::find_if( v[VertOfFace(face,i)]->GetBndVertBegin(),
                                           v[VertOfFace(face,i)]->GetBndVertEnd(),
                                           comp                                      )->GetCoord2D();

                const double f0= bdat->GetBndVal(vc2D[0]) - n_du;
                const double f1= bdat->GetBndVal(vc2D[1]) - n_du;
                const double f2= bdat->GetBndVal(vc2D[2]) - n_du;
                const double fb= bdat->GetBndVal( 1./3.*(vc2D[0] + vc2D[1] + vc2D[2]) ) - n_du;    //Barycenter of Face
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
//                const Uint v_idx=  VertOfFace(face, i);
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
//            const double absdet= FuncDet2D( v[VertOfFace(face,1)]->GetCoord() - v[VertOfFace(face,0)]->GetCoord(),
//                                            v[VertOfFace(face,2)]->GetCoord() - v[VertOfFace(face,0)]->GetCoord() );
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
    const double absdet= fabs(det);
    // Integral over tetra
    circumcircle(t, cc_center, cc_radius);
    const double P0f= Quad3CL::Quad(t, &CoeffCL::f)*6.; // absdet/vol = 6.
    _err+= 4.*cc_radius*cc_radius*P0f*P0f*absdet/6.;
//    _err+= 4.*cc_radius*cc_radius*P1DiscCL::norm_L2_sq(t, &CoeffCL::f)*absdet;

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
            const typename BndDataCL::BndSegDataCL* bdat= &Bnd.GetSegData(bidx);
            if ( bdat->IsNeumann() )
            {
//std::cerr << "neumann";
                Point2DCL vc2D[3];
                const VertexCL* v[4];

                const BndPointSegEqCL comp(bidx);
                double n_du= 0.0;  // scalar product of Du and the normal
                for (int i=0; i<4; ++i)
                {
                    v[i]= t.GetVertex(i);
                    n_du+= (v[i]->Unknowns.Exist(Idx) ? sol.Data[v[i]->Unknowns(Idx)]
                               : Bnd.GetDirBndValue(*v[i]) )
                          *(H(0,i)*normal[0] + H(1,i)*normal[1] + H(2,i)*normal[2]);
                }
                for (int i=0; i<3; ++i)
                    vc2D[i]= std::find_if( v[VertOfFace(face,i)]->GetBndVertBegin(),
                                           v[VertOfFace(face,i)]->GetBndVertEnd(),
                                           comp                                      )->GetCoord2D();

                const double f0= bdat->GetBndVal(vc2D[0]) - n_du;
                const double f1= bdat->GetBndVal(vc2D[1]) - n_du;
                const double f2= bdat->GetBndVal(vc2D[2]) - n_du;
                const double fb= bdat->GetBndVal( 1./3.*(vc2D[0] + vc2D[1] + vc2D[2]) ) - n_du;    //Barycenter of Face
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
//                const Uint v_idx=  VertOfFace(face, i);
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
//            const double absdet= FuncDet2D( v[VertOfFace(face,1)]->GetCoord() - v[VertOfFace(face,0)]->GetCoord(),
//                                            v[VertOfFace(face,2)]->GetCoord() - v[VertOfFace(face,0)]->GetCoord() );
            _err+= /* 0.5 * 2.* */cc_rad_Face * jump*jump * absdet2D/2.;
        }
    }
    return 4.*cc_radius*cc_radius*_err;
}

// PoissonP2CL
template<class Coeff>
void PoissonP2CL<Coeff>::CreateNumbering(Uint level, IdxDescCL* idx)
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
    CreateNumbOnEdge( idxnum, idx->NumUnknowns, idx->NumUnknownsEdge,
                      _MG.GetTriangEdgeBegin(level), _MG.GetTriangEdgeEnd(level), GetBndData() );
}


template<class Coeff>
void PoissonP2CL<Coeff>::DeleteNumbering(IdxDescCL* idx)
{
    const Uint idxnum = idx->GetIdx();    // idx is the index in UnknownIdxCL
    const Uint level  = idx->TriangLevel;
    idx->NumUnknowns = 0;

    // delete memory allocated for indices
    DeleteNumbOnSimplex( idxnum, _MG.GetTriangVertexBegin(level), _MG.GetTriangVertexEnd(level) );
    DeleteNumbOnSimplex( idxnum, _MG.GetTriangEdgeBegin(level), _MG.GetTriangEdgeEnd(level) );
}

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


inline double Quad( const TetraCL& t, scalar_fun_ptr f, int i, int j)
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
    double sum= a[4]*f(GetBaryCenter(t));
    for(Uint i=0; i<4; ++i)
        sum+= a[i]*f(t.GetVertex(i)->GetCoord());
    return sum;
}


inline double Quad( const TetraCL& t, scalar_fun_ptr coeff, int i)
{
    double f[5];
    
    if (i<4) // hat function on vert
    {
        f[0]= coeff( t.GetVertex(i)->GetCoord() );
        for (int k=0, l=1; k<4; ++k)
            if (k!=i) f[l++]= coeff( t.GetVertex(k)->GetCoord() );
        f[4]= coeff( GetBaryCenter(t) );
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

        double sum= vs * coeff( GetBaryCenter(t) );
        for(int k=0; k<4; ++k)
            sum+= a[k] * coeff( t.GetVertex(k)->GetCoord() );

        return sum;
    }
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


template<class Coeff>
void PoissonP2CL<Coeff>::SetupSystem(MatDescCL& Amat, VecDescCL& b) const
// Sets up the stiffness matrix and right hand side
{
    b.Clear();
    
    const Uint lvl    = Amat.GetRowLevel(),
               idx    = Amat.RowIdx->GetIdx();
    MatrixBuilderCL A( &Amat.Data, Amat.RowIdx->NumUnknowns, Amat.ColIdx->NumUnknowns); 

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupSystem: " << Amat.ColIdx->NumUnknowns 
              << " unknowns, " << std::endl;                            
	      
// fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup[10][10];
    double det, absdet;
    double tmp;

    GetGradientsOnRef(GradRef);
   
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
	
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        }

	for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }
	
        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // negative dot-product of the gradients
                coup[i][j] = QuadGrad( Grad, i, j)*absdet;
                coup[i][j]+= Quad(*sit, &_Coeff.q, i, j)*absdet;
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
                    {
                        tmp= j<4 ? _BndData.GetDirBndValue(*sit->GetVertex(j))
                                 : _BndData.GetDirBndValue(*sit->GetEdge(j-4));
                        b.Data[Numb[i]]-=          coup[j][i] * tmp;
                    }
                }
                tmp= Quad(*sit, &_Coeff.f, i)*absdet;
                b.Data[Numb[i]]+=          tmp;

                if ( i<4 ? _BndData.IsOnNeuBnd(*sit->GetVertex(i))
                         : _BndData.IsOnNeuBnd(*sit->GetEdge(i-4)) ) // vert/edge i is on Neumann boundary
                {
                    Uint face;
                    for (int f=0; f < 3; ++f)
                    {
                        face= i<4 ? FaceOfVert(i,f) : FaceOfEdge(i-4,f);
                        if ( sit->IsBndSeg(face))
                        {   
                            b.Data[Numb[i]]+=          tmp;
                        }
                    }
                }
            }      
    }
    std::cerr << "done: value part fill" << std::endl;
    
    A.Build();
}


template<class Coeff>
void PoissonP2CL<Coeff>::SetupStiffnessMatrix(MatDescCL& Amat) const
// Sets up the stiffness matrix 
{
    const Uint lvl    = Amat.GetRowLevel(),
               idx    = Amat.RowIdx->GetIdx();
    MatrixBuilderCL A( &Amat.Data, Amat.RowIdx->NumUnknowns, Amat.ColIdx->NumUnknowns); 

    IdxT Numb[10];
    bool IsOnDirBnd[10];
    
    std::cerr << "entering SetupStiffnessMatrix: " << Amat.ColIdx->NumUnknowns 
              << " unknowns, " << std::endl;                            
	      
// fill value part of matrices
    SMatrixCL<3,5> Grad[10], GradRef[10];  // jeweils Werte des Gradienten in 5 Stuetzstellen
    SMatrixCL<3,3> T;
    double coup[10][10];
    double det, absdet;
    double tmp;

    GetGradientsOnRef(GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        GetTrafoTr(T,det,*sit);
        MakeGradients(Grad, GradRef, T);
        absdet= fabs(det);
	
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and IsOnDirBnd
        for(int i=0; i<4; ++i)
        {
            if(!(IsOnDirBnd[i]= _BndData.IsOnDirBnd( *sit->GetVertex(i) )))
                Numb[i]= sit->GetVertex(i)->Unknowns(idx);
        }

	for(int i=0; i<6; ++i)
        {
            if (!(IsOnDirBnd[i+4]= _BndData.IsOnDirBnd( *sit->GetEdge(i) )))
                Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }

        // compute all couplings between HatFunctions on edges and verts
        for(int i=0; i<10; ++i)
            for(int j=0; j<=i; ++j)
            {
                // negative dot-product of the gradients
                coup[i][j] = QuadGrad( Grad, i, j)*absdet;
                coup[i][j]+= Quad(*sit, &_Coeff.q, i, j)*absdet;
                coup[j][i] = coup[i][j];
            }

        // assemble  

        for(int i=0; i<10; ++i)  // assemble row Numb[i]
            if (!IsOnDirBnd[i])  // vert/edge i is not on a Dirichlet boundary
            {
                for(int j=0; j<10; ++j)
                {
                    if (!IsOnDirBnd[j]) // vert/edge j is not on a Dirichlet boundary
                    {
                        A ( Numb[i], Numb[j] ) += coup[j][i]; 
                    }
                }

            }      
    }
    std::cerr << "done: value part fill" << std::endl;
    
    A.Build();
}


template<class Coeff>
double PoissonP2CL<Coeff>::CheckSolution(const VecDescCL& lsg, scalar_fun_ptr Lsg) const
{
    double diff, maxdiff=0, norm2= 0, L2=0;
    Uint lvl=lsg.GetLevel(),
         Idx=lsg.RowIdx->GetIdx();
    
    DiscSolCL sol(&lsg, &GetBndData(), &GetMG());
    
    std::cerr << "Abweichung von der tatsaechlichen Loesung:" << std::endl;

    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit != send; ++sit)
    {
        double absdet= sit->GetVolume()*6.,
            sum= 0;
        for(Uint i=0; i<4; ++i)
        {
            diff= (sol.val(*sit->GetVertex(i)) - Lsg(sit->GetVertex(i)->GetCoord()));
            sum+= diff*diff;
        }
        sum/= 120;
        diff= sol.val(*sit, 0.25, 0.25, 0.25) - Lsg(GetBaryCenter(*sit));
        sum+= 2./15. * diff*diff;
        L2+= sum*absdet;
    }
    L2= sqrt(L2);

    for (MultiGridCL::const_TriangVertexIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangVertexBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangVertexEnd(lvl);
         sit != send; ++sit)
    {
        if (sit->Unknowns.Exist(idx))
        {
           diff= fabs( Lsg(sit->GetCoord()) - lsg.Data[sit->Unknowns(Idx)] );
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
    return L2;
}


//============================================================================
//
//                      Marker classes for adaptive refinement
//
//============================================================================


template <class _TetraEst, class _ProblemCL>
template <class _BndData, class _VD>
void PoissonErrEstCL<_TetraEst, _ProblemCL>::Init(const P1EvalCL<double, _BndData, _VD>& sol)
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
    _InitGlobErr= sqrt(_InitGlobErr);
    _ActGlobErr= _InitGlobErr;
    if (_outp)
        *_outp << "ErrorEstimator is initialized now; initial error: " << _InitGlobErr
               << " We try to reduce by a factor of " << _RelReduction
               << "." << std::endl;
}

template <class _TetraEst, class _ProblemCL>
template <class _BndData, class _VD>
bool PoissonErrEstCL<_TetraEst, _ProblemCL>::Estimate(const P1EvalCL<double, _BndData, const _VD>& sol)
{
    const MultiGridCL& mg= _Problem.GetMG();
    const typename _ProblemCL::BndDataCL& bnd= _Problem.GetBndData();
    const Uint lvl= sol.GetLevel();
    const VecDescCL& lsg= *sol.GetSolution();
    const double exp_err= _InitGlobErr*_RelReduction/sqrt(_Meas);
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
        localerr_dist= localerr/sqrt(sit->GetVolume());
        globalerr+= localerr*localerr;
        if ( sit->IsUnrefined() )
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
    globalerr= sqrt(globalerr);
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


namespace // anonymous namespace
{
    typedef std::pair<const TetraCL*, double> Err_PairT;
    typedef std::deque<Err_PairT> Err_ContCL;

    struct AccErrCL :public std::binary_function<double, const Err_PairT, double>
    {
        double operator() (double init, const Err_PairT& ep) const
            { return init + ep.second;}
    };

    struct Err_Pair_GTCL :public std::binary_function<const Err_PairT, const Err_PairT, bool>
    {
        bool operator() (const Err_PairT& ep0, const Err_PairT& ep1)
        { return ep0.second > ep1.second; }
    };

} // end of anonymous namespace

template <class _TetraEst, class _ProblemCL>
template <class _BndData, class _VD>
void DoerflerMarkCL<_TetraEst, _ProblemCL>::Init(const P1EvalCL<double, _BndData, _VD>& sol)
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
    _InitGlobErr= sqrt(_InitGlobErr);
    _ActGlobErr= _InitGlobErr;
    if (_outp)
        *_outp << "ErrorEstimator is initialized now; initial error: " << _InitGlobErr
               << " We try to reduce by a factor of " << _RelReduction
               << " by marking the tetrahedrons with largest error-estimates, until the error marked" 
               << " is at least " << _Threshold << " of the actual global error." << std::endl;
}

template <class _TetraEst, class _ProblemCL>
template <class _BndData, class _VD>
bool DoerflerMarkCL<_TetraEst, _ProblemCL>::Estimate(const P1EvalCL<double, _BndData, const _VD>& sol)
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
    const double globalerr= sqrt(globalerr_sq);
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

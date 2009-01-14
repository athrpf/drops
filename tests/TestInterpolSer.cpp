//**************************************************************************
// File:    TestInterpolPar.cpp                                            *
// Content: testing parallel repairing of P1 and P2 functions on a         *
//          changing parallel multilevel grid                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//          Oliver Fortmeier, SC RWTH Aachen                               *
// Version: 0.1                                                            *
// Date:                                                                   *
// Begin:   01. February 2007                                              *
//**************************************************************************
/// \author Oliver Fortmeier
/// \file TestInterpolPar.cpp
/// \brief testing parallel repairing of P1 and P2 functions on a changing parallel multilevel grid

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include problem class
#include "misc/problem.h"
#include "num/fe.h"
#include "num/bndData.h"
#include "out/output.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>

using std::string;

/****************************************************************************
* G L O B A L   V A R I A B L E S                                           *
****************************************************************************/
const DROPS::Uint base_ref_x=3, base_ref_y=6, base_ref_z=3;
const DROPS::Uint refineStrategy=1; // 0 no loadbalancing, 1 adaptive repart, 2 PartKWay
const int steps=4;
const double dx=3., dy=6., dz=3.;
const char line[] ="------------------------------------------------------------";
const bool writefiles=true, printNumSimplices=false;

const int MIG=0, REF=1;         // constants for choosing the the filename in PrintMG
void PrintMG(const DROPS::MultiGridCL& mg, int type)
{
    static int REFnum=0;
    static int MIGnum=0;
    char filename[30];

    if (type==REF)
        std::sprintf(filename, "output/MG_REF_%i.mg",REFnum++);
    else if (type == MIG)
        std::sprintf(filename, "output/MG_MIG_%i.mg",MIGnum++);
    std::cout << " - Writing multigrid into: " << filename<<std::endl;

    std::ofstream file(filename);

    DROPS::DumpMGCL out(mg);
    out.put(file);
    file.close();
}

/// \brief Display number of unknowns
void DisplayNumUnknowns(__UNUSED__ const DROPS::MultiGridCL& MG, const DROPS::VecDescCL& x)
/// accumulated unknowns: accumulation of the unknowns of all processors. Some unknowns
///   are counted multiple due to the overlapping of edges and vertices
/// global unknowns: all unknowns are counted just once
{
    const DROPS::Ulint acc_num_unk = x.Data.size();
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    std::cerr << "  + Number of DOF of index "<<idx<<":  "
              <<acc_num_unk<<std::endl;
}

/// \brief Mark tetrahedra around a point for refinement
bool MarkAround(DROPS::MultiGridCL& mg, const DROPS::Point3DCL& p, double rad, int maxLevel= -1)
{
    bool mod=false;
    if (maxLevel==-1)
        maxLevel = mg.GetLastLevel()+1;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
         end= mg.GetTriangTetraEnd(); it!=end; ++it)
    {
        if ((int)it->GetLevel()<maxLevel && (GetBaryCenter(*it)-p).norm()<=rad )
        {
            it->SetRegRefMark();
            mod=true;
        }
    }
    return mod;
}

/// \brief Mark tetrahedra around a point for coarsening
bool UnMarkAround(DROPS::MultiGridCL& mg, const DROPS::Point3DCL& p, double rad)
{
    bool mod=false;
    for (DROPS::MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(),
         end= mg.GetTriangTetraEnd(); it!=end; ++it)
    {
        if ((GetBaryCenter(*it)-p).norm()<=rad )
        {
            it->SetRemoveMark();
            mod=true;
        }
    }
    return mod;
}


namespace DROPS
{
double f(const Point3DCL& p, double =0)
{
    return p.norm_sq();
}

Point3DCL f_vec(const Point3DCL& p, double =0)
{
    return p;
}


void SetFunction(VecDescCL& vec, const MultiGridCL& mg)
{
    const Uint idx=vec.RowIdx->GetIdx();
    if (vec.RowIdx->NumUnknownsVertex==1){
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
            sit!=mg.GetTriangVertexEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx))
                vec.Data[sit->Unknowns(idx)]=f(sit->GetCoord());
    }
    else if (vec.RowIdx->NumUnknownsVertex==3){
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
            sit!=mg.GetTriangVertexEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=f_vec(sit->GetCoord());
                for (int i=0; i<3; ++i)
                    vec.Data[sit->Unknowns(idx)+i]=val[i];
            }

    }
    if (vec.RowIdx->NumUnknownsEdge==1){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
            sit!=mg.GetTriangEdgeEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx))
                vec.Data[sit->Unknowns(idx)] =f(GetBaryCenter(*sit));
    }
    else if (vec.RowIdx->NumUnknownsEdge==3){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
            sit!=mg.GetTriangEdgeEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=f_vec(GetBaryCenter(*sit));
                for (int i=0; i<3; ++i)
                    vec.Data[sit->Unknowns(idx)+i] =val[i];
            }
    }
}

double CheckInterPol(const VecDescCL& vec, const MultiGridCL& mg)
{
    const Uint idx=vec.RowIdx->GetIdx();
    double vert_dist=-1, edge_dist=-1, dist=-1;

    const VertexCL* max_vert;
    const EdgeCL*   max_edge;

    if (vec.RowIdx->NumUnknownsVertex==1){
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()); sit!=mg.GetTriangVertexEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                dist =  std::fabs(vec.Data[sit->Unknowns(idx)]-f(sit->GetCoord()))
                    / std::min(vec.Data[sit->Unknowns(idx)],f(sit->GetCoord()));
                if (dist>vert_dist){
                    vert_dist = dist;
                    max_vert=&(*sit);
                }
            }
        }
    }
    else{
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()); sit!=mg.GetTriangVertexEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=MakePoint3D(vec.Data[sit->Unknowns(idx)], vec.Data[sit->Unknowns(idx)+1], vec.Data[sit->Unknowns(idx)+2]);
                dist =  (val-f_vec(sit->GetCoord())).norm()
                    / std::min(val.norm(), (f_vec(sit->GetCoord()).norm()));
                if (dist>vert_dist){
                    vert_dist = dist;
                    max_vert=&(*sit);
                }
            }
        }
    }

//     std::cerr << "Idx "<<idx<<": max difference ("<<vert_dist<<") found on Vertex: "<<max_vert->GetCoord()
//               << ": "<< vec.Data[max_vert->Unknowns(idx)]<<", instead of "<<f(max_vert->GetCoord())<<std::endl;

    if (vec.RowIdx->NumUnknownsEdge==1){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()); sit!=mg.GetTriangEdgeEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                dist=  std::fabs( vec.Data[sit->Unknowns(idx)]-f(GetBaryCenter(*sit)) )
                     / std::min( vec.Data[sit->Unknowns(idx)], f(GetBaryCenter(*sit)) );
                if (dist>edge_dist){
                    edge_dist = dist;
                    max_edge=&(*sit);
                }
            }
        }
    }
    if (vec.RowIdx->NumUnknownsEdge==3){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()); sit!=mg.GetTriangEdgeEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=MakePoint3D(vec.Data[sit->Unknowns(idx)], vec.Data[sit->Unknowns(idx)+1], vec.Data[sit->Unknowns(idx)+2]);
                dist=  ( val-f_vec(GetBaryCenter(*sit)) ).norm()
                     / std::min( val.norm(), (f_vec(GetBaryCenter(*sit))).norm() );
                if (dist>edge_dist){
                    edge_dist = dist;
                    max_edge=&(*sit);
                }
            }
        }
    }

//     if (vec.RowIdx->NumUnknownsEdge){
//         std::cerr << "Idx "<<idx<<": max difference ("<<edge_dist<<") found on Edge: "<<GetBaryCenter(*max_edge)
//                   << ": "<< vec.Data[max_edge->Unknowns(idx)]<<", instead of "<<f(GetBaryCenter(*max_edge))<<std::endl;
//     }

    dist = std::max(edge_dist,vert_dist);

    return dist;
}

typedef BndDataCL<double>    MyBoundaryCL;
typedef BndDataCL<Point3DCL> MyVecBoundaryCL;

void UpdateTriang(VecDescCL& P1_dir,
                  VecDescCL& P1_neu,
                  VecDescCL& P2_dir,
                  VecDescCL& P2_neu,
                  VecDescCL& P1_mixed,
                  VecDescCL& P2_mixed,
                  const MyBoundaryCL& dirBnd,
                  const MyBoundaryCL& neuBnd,
                  const MyVecBoundaryCL& mixedBnd,
                  MultiGridCL& mg
                 )
/** Test Repairing P1 and P2 functions with Neumann and Dirichlet boundary
    conditions.
*/
{
    // Point 1: Creating the variables
    // -----------------------------------------------------------------
    // multigrid
    const Uint LevelBeforeRef=mg.GetLastLevel();


    // P1-functions (u, v and y)
    IdxDescCL loc_u_idx(1,0,0,0), *old_u_idx=P1_dir.RowIdx,  *new_u_idx=&loc_u_idx;
    VecDescCL loc_u,              *old_u=&P1_dir,            *new_u=&loc_u;
    IdxDescCL loc_v_idx(1,0,0,0), *old_v_idx=P1_neu.RowIdx,  *new_v_idx=&loc_v_idx;
    VecDescCL loc_v,              *old_v=&P1_neu,            *new_v=&loc_v;
    IdxDescCL loc_y_idx(3,0,0,0), *old_y_idx=P1_mixed.RowIdx,*new_y_idx=&loc_y_idx;
    VecDescCL loc_y,              *old_y=&P1_mixed,          *new_y=&loc_y;

    // P2-functions (w, x and z)
    IdxDescCL loc_w_idx(1,1,0,0), *old_w_idx=P2_dir.RowIdx,  *new_w_idx=&loc_w_idx;
    VecDescCL loc_w,              *old_w=&P2_dir,            *new_w=&loc_w;
    IdxDescCL loc_x_idx(1,1,0,0), *old_x_idx=P2_neu.RowIdx,  *new_x_idx=&loc_x_idx;
    VecDescCL loc_x,              *old_x=&P2_neu,            *new_x=&loc_x;
    IdxDescCL loc_z_idx(3,3,0,0), *old_z_idx=P2_mixed.RowIdx,*new_z_idx=&loc_z_idx;
    VecDescCL loc_z,              *old_z=&P2_mixed,          *new_z=&loc_z;


    // Point 2: Refinement
    // -----------------------------------------------------------------
    mg.Refine();
    Uint LastLevel= mg.GetLastLevel();

    std::cerr <<"Level before Ref: "<<LevelBeforeRef<<" now "<<mg.GetLastLevel()<<std::endl;
    if (LevelBeforeRef>mg.GetLastLevel())
    {
        throw DROPSErrCL("MultiGrid Level has been decreased!");
    }

    std::cerr <<" Level before Ref: "<<LevelBeforeRef<<" now "<<mg.GetLastLevel()<<std::endl;
    std::cerr << " - Distribution of elements (after refinement and loadbalancing): level of mg: "<<mg.GetLastLevel()<<"\n";
    if (printNumSimplices) mg.SizeInfo(std::cout);

    // Point 5: Creation of numbering
    // -----------------------------------------------------------------
    CreateNumb(LastLevel, *new_u_idx, mg, dirBnd);
    CreateNumb(LastLevel, *new_v_idx, mg, neuBnd);
    CreateNumb(LastLevel, *new_w_idx, mg, dirBnd);
    CreateNumb(LastLevel, *new_x_idx, mg, neuBnd);
    CreateNumb(LastLevel, *new_y_idx, mg, mixedBnd);
    CreateNumb(LastLevel, *new_z_idx, mg, mixedBnd);

    // Point 6: Allocate memory for vectors
    // -----------------------------------------------------------------
    new_u->SetIdx(new_u_idx);
    new_v->SetIdx(new_v_idx);
    new_w->SetIdx(new_w_idx);
    new_x->SetIdx(new_x_idx);
    new_y->SetIdx(new_y_idx);
    new_z->SetIdx(new_z_idx);

    std::cerr <<"Number of unknowns:"<<std::endl;
    DisplayNumUnknowns(mg, *new_u);
    DisplayNumUnknowns(mg, *new_v);
    DisplayNumUnknowns(mg, *new_w);
    DisplayNumUnknowns(mg, *new_x);
    DisplayNumUnknowns(mg, *new_y);
    DisplayNumUnknowns(mg, *new_z);

    // Point 8: Repair the finite element functions
    // -----------------------------------------------------------------
    P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  funU(old_u, &dirBnd, &mg);
    RepairAfterRefineP1 (funU, *new_u);

    P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  funV(old_v, &neuBnd, &mg);
    RepairAfterRefineP1 (funV, *new_v);

    P1EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  funY(old_y, &mixedBnd, &mg);
    RepairAfterRefineP1 (funY, *new_y);

    P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  funW(old_w, &dirBnd, &mg);
    RepairAfterRefineP2 (funW, *new_w);

    P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  funX(old_x, &neuBnd, &mg);
    RepairAfterRefineP2 (funX, *new_x);

    P2EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  funZ(old_z, &mixedBnd, &mg);
    RepairAfterRefineP2 (funZ, *new_z);

    // Point 10: Delete old numbering
    // -----------------------------------------------------------------
    DeleteNumb(*old_u_idx, mg); old_u->Clear();
    DeleteNumb(*old_v_idx, mg); old_v->Clear();
    DeleteNumb(*old_w_idx, mg); old_w->Clear();
    DeleteNumb(*old_x_idx, mg); old_x->Clear();
    DeleteNumb(*old_y_idx, mg); old_y->Clear();
    DeleteNumb(*old_z_idx, mg); old_z->Clear();


    // Point 11: Put data into the old vectors
    // -----------------------------------------------------------------
    P1_dir.RowIdx->swap(*new_u_idx);  P1_dir.Clear(); P1_dir.Data=new_u->Data;
    P1_neu.RowIdx->swap(*new_v_idx);  P1_neu.Clear(); P1_neu.Data=new_v->Data;
    P2_dir.RowIdx->swap(*new_w_idx);  P2_dir.Clear(); P2_dir.Data=new_w->Data;
    P2_neu.RowIdx->swap(*new_x_idx);  P2_neu.Clear(); P2_neu.Data=new_x->Data;
    P1_mixed.RowIdx->swap(*new_y_idx);  P1_mixed.Clear(); P1_mixed.Data=new_y->Data;
    P2_mixed.RowIdx->swap(*new_z_idx);  P2_mixed.Clear(); P2_mixed.Data=new_z->Data;
}

double ZeroS( const Point3DCL&, double){ return 0; }
Point3DCL ZeroV( const Point3DCL&, double){ return Point3DCL(); }

void Strategy(MultiGridCL& mg)
{
    const Uint LastLevel=mg.GetLastLevel();

    IdxDescCL dirp1_idx(1,0,0,0);
    IdxDescCL neup1_idx(1,0,0,0);
    IdxDescCL mixedp1_idx(3,0,0,0);
    IdxDescCL dirp2_idx(1,1,0,0);
    IdxDescCL neup2_idx(1,1,0,0);
    IdxDescCL mixedp2_idx(3,3,0,0);

    const BndCondT dir[6]=
        {DirBC, DirBC, DirBC, DirBC, DirBC, DirBC};       // Dirichlet BC
    const BndSegDataCL<double>::bnd_val_fun bndfun[6]=
        {f,f, f,f, f,f};

    const BndCondT neu[6]=
        {Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC, Nat0BC};             // Neumann BC

    const BndCondT mixed[6]=
        {DirBC, DirBC, DirBC, Nat0BC, Nat0BC, Nat0BC};       // Dirichlet BC
    const BndSegDataCL<Point3DCL>::bnd_val_fun bndfunmixed[6]=
        {f_vec, f_vec, f_vec, ZeroV, ZeroV, ZeroV};

    const DROPS::BoundaryCL& bnd= mg.GetBnd();
    const DROPS::BndIdxT num_bnd= bnd.GetNumBndSeg();

    MyBoundaryCL dirBnd (num_bnd, dir, bndfun);
    MyBoundaryCL neuBnd (num_bnd, neu);
    MyVecBoundaryCL mixedBnd(num_bnd, mixed, bndfunmixed);

    CreateNumb(LastLevel, dirp1_idx, mg, dirBnd);
    CreateNumb(LastLevel, neup1_idx, mg, neuBnd);
    CreateNumb(LastLevel, dirp2_idx, mg, dirBnd);
    CreateNumb(LastLevel, neup2_idx, mg, neuBnd);
    CreateNumb(LastLevel, mixedp1_idx, mg, mixedBnd);
    CreateNumb(LastLevel, mixedp2_idx, mg, mixedBnd);

    VecDescCL dirp1; dirp1.SetIdx(&dirp1_idx);
    VecDescCL neup1; neup1.SetIdx(&neup1_idx);
    VecDescCL dirp2; dirp2.SetIdx(&dirp2_idx);
    VecDescCL neup2; neup2.SetIdx(&neup2_idx);
    VecDescCL mixedp1; mixedp1.SetIdx(&mixedp1_idx);
    VecDescCL mixedp2; mixedp2.SetIdx(&mixedp2_idx);

    for (int ref=1; ref<=steps; ++ref)
    {
        std::cerr << "------------------------------------------\n Interpolation step "
                  <<ref<<"\n------------------------------------------"<<std::endl;
        SetFunction(dirp1, mg);
        SetFunction(neup1, mg);
        SetFunction(dirp2, mg);
        SetFunction(neup2, mg);
        SetFunction(mixedp1, mg);
        SetFunction(mixedp2, mg);

        P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P1_dir(&dirp1, &dirBnd, &mg);
        P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P1_neu(&neup1, &neuBnd, &mg);
        P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P2_dir(&dirp2, &dirBnd,  &mg);
        P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P2_neu(&neup2, &neuBnd, &mg);
        P1EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  fun_P1_mixed(&mixedp1, &mixedBnd, &mg);
        P2EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  fun_P2_mixed(&mixedp2, &mixedBnd, &mg);

        Point3DCL e;
        e[0]=1.; e[1]=2.; e[2]=0.5;
        bool marked=false;
        switch (ref)
        {
            case 1:
                std::cerr << "Mark around "<<e<<std::endl;
                marked=MarkAround(mg, e, 0.5);
                break;
            case 2:
                std::cerr << "UnMark around "<<e<<std::endl;
                marked=UnMarkAround(mg, e, 0.7);
                break;
            case 3:
                std::cerr << "UnMark for ghost tetra kill"<<std::endl;
                // killedghost=UnMarkForGhostKill(mg, mg.GetLastLevel());
                // killedghost=GlobalOr(killedghost);
                // if (ProcCL::IamMaster() && killedghost)
                //    std::cerr << "A ghost tetra will be killed"<<std::endl;
                marked=true;
                break;
            default:
                if (ref%2==0)
                    marked=MarkAround(mg, Point3DCL(), 0.2, mg.GetLastLevel());
                else
                    marked=UnMarkAround(mg, Point3DCL(), 0.2);
                marked=true;
        }
        if (marked)
            std::cerr << " At least on tetra has been marked!\n";

        UpdateTriang(dirp1, neup1, dirp2, neup2, mixedp1, mixedp2,
                     dirBnd, neuBnd, mixedBnd,
                     mg);

        const double err_dir_p1= CheckInterPol(dirp1, mg),
                     err_neu_p1= CheckInterPol(neup1, mg),
                     err_dir_p2= CheckInterPol(dirp2, mg),
                     err_neu_p2= CheckInterPol(neup2, mg),
                     err_mixed_p1= CheckInterPol(neup1, mg),
                     err_mixed_p2= CheckInterPol(neup2, mg);

        std::cerr << "Maximal distance between function and interpolated P1/P2-functions:"
                  << "\n  +  P1 with Dirichlet       boundary conditions: "<<err_dir_p1
                  << "\n  +  P1 with Neumann         boundary conditions: "<<err_neu_p1
                  << "\n  +  P2 with Dirichlet       boundary conditions: "<<err_dir_p2
                  << "\n  +  P2 with Neumann         boundary conditions: "<<err_neu_p2
                  << "\n  +  P1 with mixed vectorial boundary conditions: "<<err_mixed_p1
                  << "\n  +  P2 with mixed vectorial boundary conditions: "<<err_mixed_p2
                  << std::endl;
    }
}
} // end of namespace DROPS


int main (int, char**)
{
  try
  {
    DROPS::Point3DCL e1(0.), e2(0.), e3(0.), orig(0.);
    e1[0]=dx; e2[1]=dy; e3[2]=dz;

    DROPS::MGBuilderCL * mgb= new DROPS::BrickBuilderCL(orig, e1, e2, e3, base_ref_x, base_ref_y, base_ref_z);

    DROPS::MultiGridCL mg(*mgb);

    std::cerr <<" Refine 2 times regular ... "<<std::endl;
    DROPS::MarkAll(mg);
    mg.Refine();

//     DROPS::MarkAll(mg);
//     mg.Refine();

    std::cerr <<" Refine 1 times around "<<orig<<" ... "<<std::endl;

    e1[0]=0.;  e1[1]=0.; e1[2]=0.5;
    MarkAround(mg, e1, 0.5);
    mg.Refine();

    std::cout << " - Distribution of elements: level of mg: "<<mg.GetLastLevel()<<"\n";
    mg.SizeInfo(std::cout);


    Strategy(mg);    // do all the stuff

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


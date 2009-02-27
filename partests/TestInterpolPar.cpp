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

// include parallel computing!
#include "parallel/parallel.h"
#include "parallel/parmultigrid.h"
#include "parallel/loadbal.h"
#include "parallel/partime.h"

 // include geometric computing
#include "geom/multigrid.h"
#include "geom/builder.h"

 // include in- and output
#include "out/output.h"
#include "out/ensightOut.h"

 // include problem class
#include "misc/problem.h"
#include "num/fe.h"
#include "num/bndData.h"

 // include standards
#include <iostream>
#include <iomanip>
#include <fstream>

using DROPS::ProcCL;
using std::string;

/****************************************************************************
* G L O B A L   V A R I A B L E S                                           *
****************************************************************************/
const std::string EnsDir= "ensight";
const std::string EnsCase="InterPol";
const DROPS::Uint base_ref_x=3, base_ref_y=6, base_ref_z=3;
const DROPS::Uint refineStrategy=2; // 0 no loadbalancing, 1 adaptive repart, 2 PartKWay
const int steps=4;
const double dx=3., dy=6., dz=3.;
const char line[] ="------------------------------------------------------------";
const bool writefiles=false, printNumSimplices=false, printEnsight=false;

enum TimePart{
    T_refine,
    T_migrate,
    T_repair_P1P2,
    T_createNum
};
DROPS::TimeStoreCL Times(4);

void SetDescriber()
{
    Times.SetDescriber(T_refine,      "Refinement of the multigrid");
    Times.SetDescriber(T_migrate,     "Migration of the multigrid");
    Times.SetDescriber(T_repair_P1P2, "Repairing P1 and P2 functions");
    Times.SetDescriber(T_createNum,   "Create numbering");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/// \brief Check parallel multigrid for errors
bool CheckParMultiGrid(DROPS::ParMultiGridCL& pmg)
{
    char dat[30];
    std::sprintf(dat,"output/sane%i.chk",DROPS::ProcCL::MyRank());
    std::ofstream check(dat);
    bool pmg_sane = pmg.IsSane(check),
    mg_sane  = pmg.GetMG().IsSane(check);
    check.close();
    return DROPS::ProcCL::Check(pmg_sane && mg_sane);
}

/// \brief Display number of unknowns
void DisplayNumUnknowns(const DROPS::MultiGridCL& MG, const DROPS::VecDescCL& x)
/// accumulated unknowns: accumulation of the unknowns of all processors. Some unknowns
///   are counted multiple due to the overlapping of edges and vertices
/// global unknowns: all unknowns are counted just once
{
    const DROPS::Ulint acc_num_unk = DROPS::ProcCL::GlobalSum(x.Data.size()),
                       glo_num_unk = x.RowIdx->GetGlobalNumUnknowns(MG);
    const DROPS::Uint  idx=x.RowIdx->GetIdx();

    if (DROPS::ProcCL::IamMaster())
        std::cerr << "  + Number of DOF of index "<<idx<<" (accumulated/global):  "
                  <<acc_num_unk<< "/" <<glo_num_unk<< std::endl;
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

/// \brief Mark all children of a ghost tetrahedron for removement
bool UnMarkForGhostKill (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// search for a ghost tetra and unmark all children
{
    int done=1;
    if (ProcCL::MyRank()!=0)
        ProcCL::Recv(&done, 1, ProcCL::MyRank()-1, 563738);
    else
        done=0;
    if (!done)
    {
        for (DROPS::MultiGridCL::const_TetraIterator  It(mg.GetTetrasBegin(maxLevel-1)),
            ItEnd(mg.GetTetrasEnd(maxLevel-1)); It!=ItEnd && !done; ++It)
        {
            if (It->IsGhost() && It->IsRegularlyRef()){
                for (DROPS::TetraCL::const_ChildPIterator ch(It->GetChildBegin()),
                    chEnd(It->GetChildEnd()); ch!=chEnd; ++ch)
                    (*ch)->SetRemoveMark();
                std::cerr << "Tetra "<<It->GetGID()<<" marked for ghost-kill by proc "<<ProcCL::MyRank()<<std::endl;
                done=1;
            }
        }
    }
    if (ProcCL::MyRank()<ProcCL::Size()-1)
        ProcCL::Send(&done, 1, ProcCL::MyRank()+1, 563738);

    return DROPS::ProcCL::GlobalOr(done);
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
    if (vec.RowIdx->NumUnknownsVertex()==1){
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
            sit!=mg.GetTriangVertexEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx))
                vec.Data[sit->Unknowns(idx)]=f(sit->GetCoord());
    }
    else if (vec.RowIdx->NumUnknownsVertex()==3){
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin());
            sit!=mg.GetTriangVertexEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=f_vec(sit->GetCoord());
                for (int i=0; i<3; ++i)
                    vec.Data[sit->Unknowns(idx)+i]=val[i];
            }

    }
    if (vec.RowIdx->NumUnknownsEdge()==1){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
            sit!=mg.GetTriangEdgeEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx))
                vec.Data[sit->Unknowns(idx)] =f(GetBaryCenter(*sit));
    }
    else if (vec.RowIdx->NumUnknownsEdge()==3){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin());
            sit!=mg.GetTriangEdgeEnd(); ++sit)
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=f_vec(GetBaryCenter(*sit));
                for (int i=0; i<3; ++i)
                    vec.Data[sit->Unknowns(idx)+i] =val[i];
            }
    }
}

double CheckInterPol(const VecDescCL& vec, const MultiGridCL& mg, bool checkunkset=false)
{
    const Uint idx=vec.RowIdx->GetIdx();
    double vert_dist=-1, edge_dist=-1, dist=-1;

    const VertexCL* max_vert;
    const EdgeCL*   max_edge;

    if (vec.RowIdx->NumUnknownsVertex()==1){
        for (MultiGridCL::const_TriangVertexIteratorCL sit(mg.GetTriangVertexBegin()); sit!=mg.GetTriangVertexEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                dist =  std::fabs(vec.Data[sit->Unknowns(idx)]-f(sit->GetCoord()))
                    / std::min(vec.Data[sit->Unknowns(idx)],f(sit->GetCoord()));
                if (dist>vert_dist){
                    vert_dist = dist;
                    max_vert=&(*sit);
                }
                if (!sit->Unknowns.UnkRecieved(idx) && checkunkset)
                    throw DROPSErrCL("CheckInterPol: No new unknowns on vertex!");
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
                if (!sit->Unknowns.UnkRecieved(idx) && checkunkset)
                    throw DROPSErrCL("CheckInterPol: No new unknowns on vertex!");
            }
        }
    }
//     std::cerr << "["<<ProcCL::MyRank()<<"] Idx "<<idx<<": max difference ("<<vert_dist<<") found on ("<<max_vert->GetGID()
//               <<") Vertex: "<<max_vert->GetCoord()
//               << ": "<< vec.Data[max_vert->Unknowns(idx)]<<", instead of "<<f(max_vert->GetCoord())<<std::endl;


    if (vec.RowIdx->NumUnknownsEdge()==1){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()); sit!=mg.GetTriangEdgeEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                dist=  std::fabs( vec.Data[sit->Unknowns(idx)]-f(GetBaryCenter(*sit)) )
                     / std::min( vec.Data[sit->Unknowns(idx)], f(GetBaryCenter(*sit)) );
                if (dist>edge_dist){
                    edge_dist = dist;
                    max_edge=&(*sit);
                }
                if (!sit->Unknowns.UnkRecieved(idx) && checkunkset)
                    throw DROPSErrCL("CheckInterPol: No new unknowns on edge!");
            }
        }
    }
    if (vec.RowIdx->NumUnknownsEdge()==3){
        for (MultiGridCL::const_TriangEdgeIteratorCL sit(mg.GetTriangEdgeBegin()); sit!=mg.GetTriangEdgeEnd(); ++sit){
            if (sit->Unknowns.Exist() && sit->Unknowns.Exist(idx)){
                Point3DCL val=MakePoint3D(vec.Data[sit->Unknowns(idx)], vec.Data[sit->Unknowns(idx)+1], vec.Data[sit->Unknowns(idx)+2]);
                dist=  ( val-f_vec(GetBaryCenter(*sit)) ).norm()
                     / std::min( val.norm(), (f_vec(GetBaryCenter(*sit))).norm() );
                if (dist>edge_dist){
                    edge_dist = dist;
                    max_edge=&(*sit);
                }
                if (!sit->Unknowns.UnkRecieved(idx) && checkunkset)
                    throw DROPSErrCL("CheckInterPol: No new unknowns on edge!");
            }
        }
    }

//     if (vec.RowIdx->NumUnknownsEdge){
//         std::cerr << "["<<ProcCL::MyRank()<<"] Idx "<<idx<<": max difference ("<<edge_dist<<") found on ("<<max_edge->GetGID()
//                   <<") Edge: "<<GetBaryCenter(*max_edge)
//                   << ": "<< vec.Data[max_edge->Unknowns(idx)]<<", instead of "<<f(GetBaryCenter(*max_edge))<<std::endl;
//     }

    dist = edge_dist>vert_dist ? edge_dist : vert_dist;

    return ProcCL::GlobalMax(dist);
}

typedef BndDataCL<double>    MyBoundaryCL;
typedef BndDataCL<Point3DCL> MyVecBoundaryCL;

void UpdateTriang(VecDescCL& P1_dir,
                  VecDescCL& P1_neu,
                  VecDescCL& P2_dir,
                  VecDescCL& P2_neu,
                  VecDescCL& P2_mixed,
                  const MyBoundaryCL& dirBnd,
                  const MyBoundaryCL& neuBnd,
                  const MyVecBoundaryCL& mixedBnd,
                  ParMultiGridCL& pmg,
                  LoadBalHandlerCL& lb
                 )
/** Test Repairing P1 and P2 functions with Neumann and Dirichlet boundary
    conditions.
*/
{
    // Point 1: Creating the variables
    // -----------------------------------------------------------------
    // multigrid
    MultiGridCL& mg=pmg.GetMG();
    const Uint LevelBeforeRef=mg.GetLastLevel();
    ParTimerCL timer;


    // P1-functions (u and v)
    IdxDescCL loc_u_idx( P1_FE), *old_u_idx=P1_dir.RowIdx,  *new_u_idx=&loc_u_idx;
    VecDescCL loc_u,             *old_u=&P1_dir,            *new_u=&loc_u;
    IdxDescCL loc_v_idx( P1_FE), *old_v_idx=P1_neu.RowIdx,  *new_v_idx=&loc_v_idx;
    VecDescCL loc_v,             *old_v=&P1_neu,            *new_v=&loc_v;

    // P2-functions (w, x and z)
    IdxDescCL loc_w_idx( P2_FE), *old_w_idx=P2_dir.RowIdx,  *new_w_idx=&loc_w_idx;
    VecDescCL loc_w,             *old_w=&P2_dir,            *new_w=&loc_w;
    IdxDescCL loc_x_idx( P2_FE), *old_x_idx=P2_neu.RowIdx,  *new_x_idx=&loc_x_idx;
    VecDescCL loc_x,             *old_x=&P2_neu,            *new_x=&loc_x;
    IdxDescCL loc_z_idx( vecP2_FE), *old_z_idx=P2_mixed.RowIdx,*new_z_idx=&loc_z_idx;
    VecDescCL loc_z,                *old_z=&P2_mixed,          *new_z=&loc_z;

    // And tell parallel multigrid about the unknowns
    pmg.AttachTo( &P1_dir, &dirBnd);
    pmg.AttachTo( &P1_neu, &neuBnd);
    pmg.AttachTo( &P2_dir, &dirBnd);
    pmg.AttachTo( &P2_neu, &neuBnd);
    pmg.AttachTo( &P2_mixed, &mixedBnd);

    // Point 2: Refinement
    // -----------------------------------------------------------------
    timer.Reset();
    pmg.Refine();
    timer.Stop(); Times.AddTime(T_refine, timer.GetMaxTime());
    Uint LastLevel= mg.GetLastLevel();
    if (writefiles) PrintMG(pmg, REF);

    // Point 3: Handle unknowns after refine
    // -----------------------------------------------------------------
    timer.Reset();
    pmg.HandleUnknownsAfterRefine();
    timer.Stop(); Times.AddTime(T_repair_P1P2, timer.GetMaxTime());

    // Point 4: Migration
    // -----------------------------------------------------------------
    timer.Reset();
    lb.DoMigration();
    timer.Stop(); Times.AddTime(T_migrate, timer.GetMaxTime());
    if (writefiles) PrintMG(pmg, MIG);

    if (ProcCL::IamMaster())
        std::cerr <<"["<<ProcCL::MyRank()<<"] Level before Ref: "<<LevelBeforeRef<<" now "<<mg.GetLastLevel()<<std::endl;
    if (LevelBeforeRef>mg.GetLastLevel())
    {
        throw DROPSErrCL("MultiGrid Level has been decreased!");
    }
    if (ProcCL::IamMaster() && printNumSimplices)
        std::cerr << " - Distribution of elements (after refinement and loadbalancing): level of mg: "<<mg.GetLastLevel()<<"\n";
    if (printNumSimplices) mg.SizeInfo(std::cout);

    // Point 5: Creation of numbering
    // -----------------------------------------------------------------
    timer.Reset();
    new_u_idx->CreateNumbering( LastLevel, mg, dirBnd);
    new_v_idx->CreateNumbering( LastLevel, mg, neuBnd);
    new_w_idx->CreateNumbering( LastLevel, mg, dirBnd);
    new_x_idx->CreateNumbering( LastLevel, mg, neuBnd);
    new_z_idx->CreateNumbering( LastLevel, mg, mixedBnd);
    timer.Stop(); Times.AddTime(T_createNum, timer.GetMaxTime());

    // Point 6: Allocate memory for vectors
    // -----------------------------------------------------------------
    new_u->SetIdx(new_u_idx);
    new_v->SetIdx(new_v_idx);
    new_w->SetIdx(new_w_idx);
    new_x->SetIdx(new_x_idx);
    new_z->SetIdx(new_z_idx);

    if (ProcCL::IamMaster())
        std::cerr <<"Number of unknowns:"<<std::endl;
    DisplayNumUnknowns(mg, *new_u);
    DisplayNumUnknowns(mg, *new_v);
    DisplayNumUnknowns(mg, *new_w);
    DisplayNumUnknowns(mg, *new_x);
    DisplayNumUnknowns(mg, *new_z);

    // Point 7: Handle transfered unknowns
    // -----------------------------------------------------------------
    timer.Reset();
    pmg.HandleNewIdx(old_u_idx, new_u);
    pmg.HandleNewIdx(old_v_idx, new_v);
    pmg.HandleNewIdx(old_w_idx, new_w);
    pmg.HandleNewIdx(old_x_idx, new_x);
    pmg.HandleNewIdx(old_z_idx, new_z);
    pmg.DeleteRecvBuffer();

    // Point 8: Repair the finite element functions
    // -----------------------------------------------------------------
    P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  funU(old_u, &dirBnd, &mg);
    RepairAfterRefineP1 (funU, *new_u);
    pmg.CompleteRepair(new_u);
    CheckInterPol(*new_u, mg, true);

    P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  funV(old_v, &neuBnd, &mg);
    RepairAfterRefineP1 (funV, *new_v);
    pmg.CompleteRepair(new_v);
    CheckInterPol(*new_u, mg, true);

    P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  funW(old_w, &dirBnd, &mg);
    RepairAfterRefineP2 (funW, *new_w);
    pmg.CompleteRepair(new_w);
    CheckInterPol(*new_u, mg, true);

    P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  funX(old_x, &neuBnd, &mg);
    RepairAfterRefineP2 (funX, *new_x);
    pmg.CompleteRepair(new_x);
    CheckInterPol(*new_u, mg, true);

    P2EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  funZ(old_z, &mixedBnd, &mg);
    RepairAfterRefineP2 (funZ, *new_z);
    pmg.CompleteRepair(new_z);
    CheckInterPol(*new_u, mg, true);

    // Point 9: Delete all flags for recieved unknowns
    // -----------------------------------------------------------------
    pmg.DelAllUnkRecv();
    pmg.DeleteRecvBuffer();
    pmg.DeleteVecDesc();

    timer.Stop(); Times.AddTime(T_repair_P1P2, timer.GetMaxTime());
    CheckParMultiGrid(pmg);

    // Point 10: Delete old numbering
    // -----------------------------------------------------------------
    old_u_idx->DeleteNumbering( mg); old_u->Clear();
    old_v_idx->DeleteNumbering( mg); old_v->Clear();
    old_w_idx->DeleteNumbering( mg); old_w->Clear();
    old_x_idx->DeleteNumbering( mg); old_w->Clear();
    old_z_idx->DeleteNumbering( mg); old_z->Clear();


    // Point 11: Put data into the old vectors
    // -----------------------------------------------------------------
    P1_dir.RowIdx->swap(*new_u_idx);  P1_dir.Clear(); P1_dir.Data=new_u->Data;
    P1_neu.RowIdx->swap(*new_v_idx);  P1_neu.Clear(); P1_neu.Data=new_v->Data;
    P2_dir.RowIdx->swap(*new_w_idx);  P2_dir.Clear(); P2_dir.Data=new_w->Data;
    P2_neu.RowIdx->swap(*new_x_idx);  P2_neu.Clear(); P2_neu.Data=new_x->Data;
    P2_mixed.RowIdx->swap(*new_z_idx);  P2_mixed.Clear(); P2_mixed.Data=new_z->Data;
}

double ZeroS( const Point3DCL&, double){ return 0; }
Point3DCL ZeroV( const Point3DCL&, double){ return Point3DCL(); }

void Strategy(ParMultiGridCL& pmg, LoadBalHandlerCL& lb)
{
    MultiGridCL& mg= pmg.GetMG();
    const Uint LastLevel=mg.GetLastLevel();

    IdxDescCL dirp1_idx(   P1_FE);
    IdxDescCL neup1_idx(   P1_FE);
    IdxDescCL dirp2_idx(   P2_FE);
    IdxDescCL neup2_idx(   P2_FE);
    IdxDescCL mixedp2_idx( vecP2_FE);

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

    dirp1_idx.CreateNumbering( LastLevel, mg, dirBnd);
    neup1_idx.CreateNumbering( LastLevel, mg, neuBnd);
    dirp2_idx.CreateNumbering( LastLevel, mg, dirBnd);
    neup2_idx.CreateNumbering( LastLevel, mg, neuBnd);
    mixedp2_idx.CreateNumbering( LastLevel, mg, mixedBnd);

    VecDescCL dirp1; dirp1.SetIdx(&dirp1_idx);
    VecDescCL neup1; neup1.SetIdx(&neup1_idx);
    VecDescCL dirp2; dirp2.SetIdx(&dirp2_idx);
    VecDescCL neup2; neup2.SetIdx(&neup2_idx);
    VecDescCL mixedp2; mixedp2.SetIdx(&mixedp2_idx);

    pmg.AttachTo(&dirp1, &dirBnd);
    pmg.AttachTo(&neup1, &neuBnd);
    pmg.AttachTo(&dirp2, &dirBnd);
    pmg.AttachTo(&neup2, &neuBnd);
    pmg.AttachTo(&mixedp2, &mixedBnd);


    // Erzeuge ensight case File und geom-File

    EnsightP2SolOutCL *ensight = new EnsightP2SolOutCL( mg, &neup2_idx, false);
    ensight->SetMasterOut();
    const string filename= EnsDir + "/" + EnsCase;
    const string datgeo=   filename+".geo";

    string *dattmp1 = new string(filename+".dat1");
    string *dattmp2 = new string(filename+".dat2");
    string *dattmp3 = new string(filename+".dat3");
    string *dattmp4 = new string(filename+".dat4");

    ensight->CaseBegin( string(EnsCase+".case").c_str(), 6);
    ensight->DescribeGeom(   "Wuerfel", datgeo, true);
    ensight->DescribeScalar( "DirP1", *dattmp1, true);
    ensight->DescribeScalar( "NeuP1", *dattmp2, true);
    ensight->DescribeScalar( "DirP2", *dattmp3, true);
    ensight->DescribeScalar( "NeuP2", *dattmp4, true);
    ensight->Commit();

    for (int ref=1; ref<=steps; ++ref)
    {
        ParTimerCL timer;
        timer.Start();
        if (ProcCL::IamMaster())
            std::cerr << "------------------------------------------\n Interpolation step "
                      <<ref<<"\n------------------------------------------"<<std::endl;
        SetFunction(dirp1, mg);
        SetFunction(neup1, mg);
        SetFunction(dirp2, mg);
        SetFunction(neup2, mg);
        SetFunction(mixedp2, mg);

        P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P1_dir(&dirp1, &dirBnd, &mg);
        P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P1_neu(&neup1, &neuBnd, &mg);
        P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P2_dir(&dirp2, &dirBnd,  &mg);
        P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P2_neu(&neup2, &neuBnd, &mg);
        P2EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  fun_P2_mixed(&mixedp2, &mixedBnd, &mg);

        if (printEnsight){
            if (ProcCL::IamMaster())
                std::cerr << "Write ensightfiles for time "<<ref<<std::endl;

            ensight->putGeom( datgeo, ref);
            ensight->putScalar(*dattmp1, fun_P1_dir, ref);
            ensight->putScalar(*dattmp2, fun_P1_neu, ref);
            ensight->putScalar(*dattmp3, fun_P2_dir, ref);
            ensight->putScalar(*dattmp4, fun_P2_neu, ref);
            ensight->Commit();
        }

        Point3DCL e;
        e[0]=1.; e[1]=2.; e[2]=0.5;
        bool marked=false;
        switch (ref)
        {
            case 1:
                if (ProcCL::IamMaster())
                    std::cerr << "Mark around "<<e<<std::endl;
                marked=MarkAround(mg, e, 0.5);
                break;
            case 2:
                if (ProcCL::IamMaster())
                    std::cerr << "UnMark around "<<e<<std::endl;
                marked=UnMarkAround(mg, e, 0.7);
                break;
            case 3:
//                 if (ProcCL::IamMaster())
//                     std::cerr << "UnMark for ghost tetra kill"<<std::endl;
//                 killedghost=UnMarkForGhostKill(mg, mg.GetLastLevel());
//                 killedghost=GlobalOr(killedghost);
//                 marked=killedghost;
//                 if (ProcCL::IamMaster() && killedghost)
//                    std::cerr << "A ghost tetra will be killed"<<std::endl;
                marked=true;
                break;
            default:
                if (ref%2==0)
                    marked=MarkAround(mg, Point3DCL(), 0.2, mg.GetLastLevel());
                else
                    marked=UnMarkAround(mg, Point3DCL(), 0.2);
                marked=true;
        }
        marked= ProcCL::GlobalOr(marked);
        if (ProcCL::IamMaster() && marked)
            std::cerr << " At least on tetra has been marked!\n";

        UpdateTriang(dirp1, neup1, dirp2, neup2, mixedp2,
                     dirBnd, neuBnd, mixedBnd,
                     pmg, lb);

        const double err_dir_p1= CheckInterPol(dirp1, mg),
                     err_neu_p1= CheckInterPol(neup1, mg),
                     err_dir_p2= CheckInterPol(dirp2, mg),
                     err_neu_p2= CheckInterPol(neup2, mg),
                     err_mixed_p2= CheckInterPol(neup2, mg);

        if (ProcCL::IamMaster())
            std::cerr << "Maximal distance between function and interpolated P1/P2-functions:"
                      << "\n  +  P1 with Dirichlet       boundary conditions: "<<err_dir_p1
                      << "\n  +  P1 with Neumann         boundary conditions: "<<err_neu_p1
                      << "\n  +  P2 with Dirichlet       boundary conditions: "<<err_dir_p2
                      << "\n  +  P2 with Neumann         boundary conditions: "<<err_neu_p2
                      << "\n  +  P2 with mixed vectorial boundary conditions: "<<err_mixed_p2
                      << std::endl;

        P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P1_dir2(&dirp1, &dirBnd, &mg);
        P1EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P1_neu2(&neup1, &neuBnd, &mg);
        P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P2_dir2(&dirp2, &dirBnd,  &mg);
        P2EvalCL<double, const MyBoundaryCL, const VecDescCL>  fun_P2_neu2(&neup2, &neuBnd, &mg);
        P2EvalCL<Point3DCL, const MyVecBoundaryCL, const VecDescCL>  fun_P2_mixed2(&mixedp2, &mixedBnd, &mg);


        if (printEnsight){
            if (ProcCL::IamMaster())
                std::cerr << "Write ensightfiles for time "<<(0.5+ref)<<std::endl;

            ensight->putGeom( datgeo, 0.5+ref);
            ensight->putScalar(*dattmp1, fun_P1_dir2, 0.5+ref);
            ensight->putScalar(*dattmp2, fun_P1_neu2, 0.5+ref);
            ensight->putScalar(*dattmp3, fun_P2_dir2, 0.5+ref);
            ensight->putScalar(*dattmp4, fun_P2_neu2, 0.5+ref);
            ensight->Commit();
        }

        timer.Stop();
        double duration= timer.GetMaxTime();
        if (ProcCL::IamMaster())
            std::cerr << "Step "<<ref<<" took "<<duration<<" sec"<<std::endl;
    }

    ensight->CaseEnd();
}
} // end of namespace DROPS


int main (int argc, char** argv)
{
  DROPS::ProcInitCL procinit(&argc, &argv);
  DROPS::ParMultiGridInitCL pmginit;
  try
  {

    DROPS::ParTimerCL alltime;
    SetDescriber();
    DROPS::ParMultiGridCL pmg= DROPS::ParMultiGridCL::Instance();

    DROPS::Point3DCL e1(0.), e2(0.), e3(0.), orig(0.);
    e1[0]=dx; e2[1]=dy; e3[2]=dz;

    DROPS::ParTimerCL time;
    DROPS::MGBuilderCL * mgb;
    if (DROPS::ProcCL::IamMaster())
        mgb = new DROPS::BrickBuilderCL(orig, e1, e2, e3, base_ref_x, base_ref_y, base_ref_z);
    else
        mgb = new DROPS::EmptyBrickBuilderCL(orig, e1, e2, e3);


    DROPS::MultiGridCL mg(*mgb);
    pmg.AttachTo(mg);

    DROPS::LoadBalHandlerCL lb(mg);
    lb.DoInitDistribution(DROPS::ProcCL::Master());
    switch (refineStrategy){
        case 0 : lb.SetStrategy(DROPS::NoMig);     break;
        case 1 : lb.SetStrategy(DROPS::Adaptive);  break;
        case 2 : lb.SetStrategy(DROPS::Recursive); break;
    }

    if (DROPS::ProcCL::IamMaster())
        std::cerr <<" Refine 2 times regular ... "<<std::endl;
    DROPS::MarkAll(mg);
    pmg.Refine();
    lb.DoMigration();

    DROPS::MarkAll(mg);
    pmg.Refine();
    lb.DoMigration();


    if (DROPS::ProcCL::IamMaster())
        std::cerr <<" Refine 1 times around "<<orig<<" ... "<<std::endl;

    e1[0]=0.;  e1[1]=0.; e1[2]=0.5;
    MarkAround(mg, e1, 0.5);
    pmg.Refine();
    if (writefiles) DROPS::PrintMG(pmg, DROPS::REF);
    lb.DoMigration();
    if (writefiles) DROPS::PrintMG(pmg, DROPS::MIG);

    MarkAround(mg, e1, 0.5);
    pmg.Refine();
    lb.DoMigration();

    if (DROPS::ProcCL::IamMaster())
        std::cout << " - Distribution of elements: level of mg: "<<mg.GetLastLevel()<<"\n";
    mg.SizeInfo(std::cout);


    Strategy(pmg, lb);    // do all the stuff

    alltime.Stop();
    Times.SetOverall(alltime.GetMaxTime());
    Times.Print(std::cout);

    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}


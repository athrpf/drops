/// \file TestRefPar.cpp
/// \brief Testing parallel refinement and coarsening
/// \author LNM RWTH Aachen: ; SC RWTH Aachen: Oliver Fortmeier

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

// "parallel" Header-Files
#include "parallel/parmultigrid.h"
#include "parallel/partime.h"
#include "parallel/metispartioner.h"
#include "parallel/loadbal.h"
#include "parallel/parmgserialization.h"

// geometry Header-Files
#include "geom/builder.h"
#include "geom/multigrid.h"

// Ausgabe in geomview-format
#include "out/output.h"
#include "partests/params.h"

// Standard-Header-Files fuer Ausgaben
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

using namespace std;

/****************************************************************************
* G L O B A L E  V A R I A B L E N                                          *
****************************************************************************/
// Zeiten, die gemessen werden sollen
enum TimePart{
    T_Ref,
    T_SetupGraph,
    T_CalcDist,
    T_Migration,
    T_Check,
    T_print
};

// Tabelle, in der die Zeiten stehen
DROPS::TimeStoreCL Times(6);

// Parameter-Klasse
DROPS::ParamParRefCL C;

/****************************************************************************
    * S E T   D E S C R I B E R   F O R   T I M E S T O R E  C L                *
****************************************************************************/
void SetDescriber()
{
    Times.SetDescriber(T_Ref, "Refinement");
    Times.SetDescriber(T_SetupGraph, "Setup loadbalancing graph");
    Times.SetDescriber(T_CalcDist, "Calculate Distribution");
    Times.SetDescriber(T_Migration, "Migration");
    Times.SetDescriber(T_Check, "Checking MG");
    Times.SetDescriber(T_print, "Printing");
    Times.SetCounterDescriber("Moved MultiNodes");
}

/****************************************************************************
* F I L E  H A N D L I N G                                                  *
****************************************************************************/
void PrintGEO(const DROPS::ParMultiGridCL& pmg)
{
    const int me=DROPS::ProcCL::MyRank();
    static int num=0;
    char filename[30];
    sprintf (filename, "output/geo_%i_GEOM_%i.geo",me,num++);
    ofstream file(filename);
    file << DROPS::GeomMGOutCL(pmg.GetMG(),-1,false,0.08,0.15) << std::endl;
    file.close();
}

/****************************************************************************
* C H E C K  P A R  M U L T I  G R I D                                      *
*****************************************************************************
*   Checkt, ob die parallele Verteilung und die MultiGrid-Struktur gesund   *
*   zu sein scheint. Zudem wird der Check von DDD aufgerufen                *
****************************************************************************/
void CheckParMultiGrid(DROPS::ParMultiGridCL& pmg, int type, int proc=0)
{
    const int me=DROPS::ProcCL::MyRank();
    DROPS::ParTimerCL time;
    double duration;

    if (type==DROPS::MIG && !C.checkMig) return;
    if (type==DROPS::REF && !C.checkRef) return;
    if (me==proc)
        std::cout << "  - Check des parallelen MultiGrids ... ";

    char dat[30];
    sprintf(dat,"output/sane%i.chk",me);
    ofstream check(dat);
    time.Reset();
    bool pmg_sane = pmg.IsSane(check),
         mg_sane  = pmg.GetMG().IsSane(check);

    if (DROPS::ProcCL::Check(pmg_sane && mg_sane)){
         if (me==proc) std::cout << "OK\n";
    }
    else{
        // Always exit on error
        if (me==proc) std::cout << " nicht OK!!!\n";
            std::cout << "EXIT Error found in multigrid\n";
            exit(-1);
    }
    if (C.checkDDD){
        if (me==proc) cout << "  -";
        pmg.ConsCheck();
    }
    time.Stop();
    duration = time.GetMaxTime();
    Times.AddTime(T_Check,duration);
    if (C.printTime && me==proc) std::cout << "       --> "<<duration<<" sec\n";
    check.close();
}

/****************************************************************************
* M A R K I N G   S T R A T E G I E S                                       *
*****************************************************************************
*   Setze Markierungen auf den Tetraedern. Zum Verfeinern und Vergroebern.  *
****************************************************************************/
void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}

void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

void MarkCorner (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Corner; Corner[0]=0.; Corner[1]=0.; Corner[2]=0.;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
         ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Corner).norm()<=0.3)
            It->SetRegRefMark();
    }
}

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

bool UnMarkForGhostKill (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
// search for a ghost tetra and unmark all children
{
    int done=1;
    if (DROPS::ProcCL::MyRank()!=0)
        DROPS::ProcCL::Recv(&done, 1, DROPS::ProcCL::MyRank()-1, 563738);
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
                std::cout << "Tetra "<<It->GetGID()<<" marked for ghost-kill by proc "<<DROPS::ProcCL::MyRank()<<std::endl;
                done=1;
            }
        }
    }
    if (DROPS::ProcCL::MyRank()<DROPS::ProcCL::Size()-1)
        DROPS::ProcCL::Send(&done, 1, DROPS::ProcCL::MyRank()+1, 563738);


    return DROPS::ProcCL::GlobalOr(done);
}

void CreateInitGrid(DROPS::ParMultiGridCL& pmg,  int proc)
{
    using namespace DROPS;
    MultiGridCL *mg;
    const int me = ProcCL::MyRank(), size = ProcCL::Size();
    DROPS::ParTimerCL time;
    double duration;

    Point3DCL e1(0.0), e2(0.0), e3(0.0), orig(0.0);

    if(C.init_cond==0)
    {
        e1[0]=C.brk_dim[0]; e2[1]=C.brk_dim[1]; e3[2]= C.brk_dim[2];
        if (ProcCL::MyRank()==proc)
        {
            BrickBuilderCL brick(C.brk_orig, e1, e2, e3, C.brk_BasicRefX, C.brk_BasicRefY, C.brk_BasicRefZ);
            mg = new DROPS::MultiGridCL(brick);
        }
        else
        {
            EmptyBrickBuilderCL emptyBrick(C.brk_orig, e1, e2, e3, C.brk_BasicRefX);
            mg = new DROPS::MultiGridCL(emptyBrick);
        }
        pmg.AttachTo(*mg);

        LoadBalCL lb(*mg);
        if (size>1)
        {
            lb.DeleteGraph();
            if (me==0)
                std::cout << "  - Erstelle dualen reduzierten Graphen ...\n";
            time.Reset();
            lb.CreateDualRedGraph(true);
            time.Stop();
            if (C.printTime){
                duration = time.GetMaxTime();
                if (me==0) std::cout << "       --> "<<duration<<" sec\n";
            }

            if (me==0)
                std::cout << "  - Berechne eine Graphpartitionierung ...\n";
            time.Reset();
            lb.PartitionSer(proc);
            time.Stop();
            if (C.printTime){
                duration = time.GetMaxTime();
                if (me==0) std::cout << "       --> "<<duration<<" sec\n";
            }

            if (me==0)
                std::cout << "  - Migration ...\n";

            // Verteile das MultiGrid
            time.Reset();
            pmg.XferStart();
            lb.Migrate();
            pmg.XferEnd();
            time.Stop();
            if (C.printTime){
                duration = time.GetMaxTime();
                if (me==0) std::cout << "       --> "<<duration<<" sec\n";
            }
            Times.IncCounter(lb.GetMovedMultiNodes());
            //pmg.DelAllUnkRecv();
        }
    }

    else if (C.init_cond==1)
    {
        e1[0]=e2[1]=e3[2]= 1.;
        if (ProcCL::MyRank()==proc)
        {
            BrickBuilderCL builder(orig, e1, e2, e3, 4, 4, 4);
            FileBuilderCL fileBuilder(C.init_pre, &builder);
            mg = new DROPS::MultiGridCL(fileBuilder);

            std::ofstream serSanity("sanity.txt");

            std::cout << "\n \n MultiGrid mit "<<mg->GetNumLevel()<<" Leveln aus Datei gelesen\n \n";
            serSanity << SanityMGOutCL(*mg) << std::endl;
        }
        else
        {
            EmptyBrickBuilderCL builder(orig, e1, e2, e3, C.refined+1);
            mg = new DROPS::MultiGridCL(builder);
        }

        pmg.AttachTo(*mg);

        LoadBalHandlerCL lb(*mg);
        lb.DoInitDistribution();
    }

    else
    {
        throw DROPSErrCL("Unknown init condition");
    }
}

void DoMigration(DROPS::ParMultiGridCL &pmg, DROPS::LoadBalCL &LoadBal, int lb)
{
    DROPS::ParTimerCL time;
    double duration;
    const int me = DROPS::ProcCL::MyRank(), size = DROPS::ProcCL::Size();
    if (size>0 && lb!=0)
    {
        if (me==0) cout << "  - Erstelle Graphen ... \n";
        LoadBal.DeleteGraph();
        time.Reset();
        LoadBal.CreateDualRedGraph();
        time.Stop(); duration = time.GetMaxTime();
        Times.AddTime(T_SetupGraph, duration);
        if (C.printTime && me==0) std::cout << "       --> "<<duration<<" sec\n";

        if (me==0) cout << "  - Erstelle Partitionen ... \n";
        time.Reset();
        LoadBal.PartitionPar();

        time.Stop(); duration = time.GetMaxTime();
        Times.AddTime(T_CalcDist, duration);
        if (C.printTime && me==0) std::cout << "       --> "<<duration<<" sec\n";
        if (me==0) cout << "  - Migration ... \n";
        time.Reset();
        if (lb!=0){
            pmg.XferStart();
            LoadBal.Migrate();
            pmg.XferEnd();
            //pmg.DelAllUnkRecv();
        }
        time.Stop(); duration = time.GetMaxTime();
        Times.AddTime(T_Migration, duration);
        if (C.printTime && me==0)
            std::cout << "       --> "<<duration<<" sec, moved MulitNodes " << LoadBal.GetMovedMultiNodes() << "\n";
        Times.IncCounter(LoadBal.GetMovedMultiNodes());
    }
}

using namespace DROPS;
/****************************************************************************
* M A I N                                                                   *
****************************************************************************/
int main(int argc, char* argv[])
{
    DROPS::ProcInitCL procinit(&argc, &argv);
    DROPS::ParMultiGridInitCL pmginit;
    if (C.printTime)
        DROPS::ParTimerCL::TestBandwidth(std::cout);
    try
    {
        const char line[] = "----------------------------------------------------------------------------------";
        const char dline[]= "==================================================================================";
        SetDescriber();

        //DDD_SetOption(OPT_XFER_PRUNE_DELETE,OPT_ON);

        DROPS::ParTimerCL alltime, time;
        double duration;

        const int me     = DROPS::ProcCL::MyRank();
        const int master = 0;

        // Parameter file einlesen ...
        if (argc!=2){
            std::cout << "You have to specify one parameter:\n\t" << argv[0] << " <param_file>" << std::endl; return 1;
        }
        std::ifstream param( argv[1]);
        if (!param){
            std::cout << "error while opening parameter file\n"; return 1;
        }

        param >> C;

        param.close();
        if (me==master)
            std::cout << C << std::endl;

        if (me==0){
            cout << dline << endl << " + Erstelle initiale Grid (Wuerfel der Laenge 1) auf Prozessor " << master <<": ...\n";
        }
        DROPS::ParMultiGridCL pmg = DROPS::ParMultiGridCL::Instance();
        pmg.ShowTypes();
        //pmg.ShowInterfaces();
        CreateInitGrid(pmg, master);

        DROPS::MultiGridCL &mg = pmg.GetMG();
        DROPS::LoadBalCL LoadBal(mg);

        if (C.printSize){
            if (me==master) cout << "  - Verteilung der Elemente:\n";
            mg.SizeInfo(cout);
        }
        if (C.printPMG){
            if (me==master) cout << " + Schreibe Debug-Informationen in ein File ... ";
            PrintMG(pmg);
            if (me==master)  cout << " OK\n";
        }
        if (C.printGEO){
            if (me==master) cout << " + Schreibe das Multigrid im Geomview-Format in ein File ... ";
            PrintGEO(pmg);
            if (me==master) cout << " OK\n";
        }

        CheckParMultiGrid(pmg,REF,master);

        if (me==master){
            cout << dline << endl << " Verfeinere das Gitter nun " << C.markall << " mal global, " << C.markdrop
                    << " mal in der Mitte um den Tropfen\n und " << C.markcorner << " mal um der Ecke (0,0,0)\n"
                    << " Es wird die Strategie ";
            switch (C.refineStrategy){
                case 0 : cout << "No Loadbalancing ";break;
                case 1 : cout << "AdaptiveRefine "; break;
                case 2 : cout << "PartKWay "; break;
                default: cout << "Unbekannte Strategy ...\n EXIT"; exit(0);
            }
            cout << "verwendet. Es markiert der Prozessor " << C.markingproc << "\n" << dline << endl;
        }
        int movedRefNodes=0, movedCoarseNodes=0;
        int numrefs;

        switch (C.Strategy)
        {
            case 0:  numrefs= C.markall+C.markdrop+C.markcorner; break;
            case 1:  numrefs=5; break;
            case 2:  numrefs=C.markall; break;
            case 3:  numrefs=C.markall; break;
            default: throw DROPSErrCL("Specify the refinement strategy!");
        }

        for (int ref=0; ref<C.markall+C.markdrop+C.markcorner; ++ref)
        {
            DROPS::Point3DCL e, e1;
            bool marked=false;
            bool killedghost=false;

            switch (C.Strategy)
            {
            case 0:
                if (me==master) cout << " + Refine " << (ref) << " : ";
                if (ref < C.markall){
                    if (me==0) cout << "all ...\n";
                    if (C.markingproc==-1 || C.markingproc==me)
                        DROPS::MarkAll(mg);
                }
                else if (ref < C.markdrop+C.markall){
                    if (me==master) cout << "drop ...\n";
                    if (C.markingproc==-1 || C.markingproc==me)
                        MarkDrop(mg, mg.GetLastLevel());
                }
                else{
                    if (me==master) cout << "corner ...\n";
                    if (C.markingproc==-1 || C.markingproc==me)
                        MarkCorner(mg, mg.GetLastLevel());
                }
            break;
            case 1:
                e[0]=1.; e[1]=2.; e[2]=0.5;
                e1[0]=0.;  e1[1]=0.; e1[2]=0.5;
                switch (ref)
                {
                    case 0:
                        if (ProcCL::IamMaster())
                            std::cout << "Mark all "<<std::endl;
                        MarkAll(mg);
                        marked=true;
                        break;
                    case 1:
                        if (ProcCL::IamMaster())
                            std::cout << "Mark all"<<std::endl;
                        MarkAll(mg);
    //                     marked=MarkAround(mg, e1, 0.5);
                        marked=true;
                        break;
                    case 2:
                        if (ProcCL::IamMaster())
                            std::cout << "Mark around "<<e<<std::endl;
                        marked=MarkAround(mg, e, 0.5);
                        break;
                    case 3:
                        if (ProcCL::IamMaster())
                            std::cout << "UnMark around "<<e<<std::endl;
                        marked=UnMarkAround(mg, e, 0.6);
                        break;
                    case 4:
                        if (ProcCL::IamMaster())
                            std::cout << "UnMark for ghost tetra kill"<<std::endl;
                        killedghost=UnMarkForGhostKill(mg, mg.GetLastLevel());
                        killedghost= ProcCL::GlobalOr(killedghost);
                        if (ProcCL::IamMaster() && killedghost)
                            std::cout << "A ghost tetra will be killed"<<std::endl;
                        break;
                    default:
                        std::cout << "I do not know this case!\n";
                }
            break;
            case 3:
                if (ref%2==0)
                    MarkAll(mg);
                else
                    UnMarkAll(mg);
                marked=true;
            break;
            }   // end of switch C.Strategy


            time.Reset(); pmg.Refine(); time.Stop();
            duration = time.GetMaxTime();
            Times.AddTime(T_Ref,duration);
            if (C.printTime && me==master) std::cout << "       --> "<<duration<<" sec\n";

            if (C.printPMG){
                if (me==master) cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg,REF);
                if (me==master) cout << " OK\n";
            }

            if (C.checkRef)
                CheckParMultiGrid(pmg,REF,master);

            DynamicDataInterfaceCL::ConsCheck();
            DoMigration(pmg, LoadBal,C.refineStrategy);
            movedRefNodes += LoadBal.GetMovedMultiNodes();

            if (C.printPMG){
                if (me==master) cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg,MIG);
                if (me==master) cout << " OK\n";
            }
            if (C.printGEO){
                if (me==master) cout << "  - Schreibe das Multigrid im Geomview-Format in ein File ... ";
                PrintGEO(pmg);
                if (me==master) cout << " OK\n";
            }
            if (C.printSize){
                if (me==master)
                    cout << "  - Verteilung der Elemente:\n";
                mg.SizeInfo(cout);
            }

            if (C.checkMig)
                CheckParMultiGrid(pmg,MIG,master);

            if (me==master && ref!=C.markall+C.markdrop+C.markcorner-1) cout << line << endl;
        }

        if (C.middleMig){
            if (me==master)
                cout <<dline<<endl<< " + Last-Verteilung zwischen dem Verfeinern und Vergroebern ...\n";
            DoMigration(pmg, LoadBal,0);
            movedRefNodes += LoadBal.GetMovedMultiNodes();
            CheckParMultiGrid(pmg,MIG,master);
        }

        if (me==master){
            cout <<dline<<endl << " Vergroebere nun das Gitter zunaechst " << C.coarsedrop
                    << " mal um den Tropfen herum und dann " << C.coarseall << " ueberall\n Es wird die Strategie ";
            switch (C.coarseStrategy){
                case 0 : cout << "No Loadbalancing ";break;
                case 1 : cout << "AdaptiveRefine "; break;
                case 2 : cout << "PartKWay "; break;
                default: cout << "Unbekannte Strategy ...\n EXIT"; exit(0);
            }
            cout << "verwendet. Es markiert der Prozessor " << C.unmarkingproc << "\n" << dline << endl;
        }

        for (int ref =0; ref<C.coarsedrop+C.coarseall; ++ref)
        {
            if (ref < C.coarsedrop){
                cout << " + Coarse drop (" << ref << ") ... \n";
                if (C.unmarkingproc==-1 || C.unmarkingproc==me)
                    UnMarkDrop(mg, mg.GetLastLevel());
            }
            else {
                cout << " + Coarse all (" << ref << ") ... \n";
                if (C.unmarkingproc==-1 || C.unmarkingproc==me){
                    DROPS::UnMarkAll(mg);
                }
            }

            time.Reset(); pmg.Refine(); time.Stop();
            duration = time.GetMaxTime();
            Times.AddTime(T_Ref,duration);
            if (C.printTime && me==master) std::cout << "       --> "<<duration<<" sec\n";
            if (C.printPMG)
            {
                if (me==master) cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg, REF);
                if (me==master) cout << " OK\n";
            }
            if (C.printGEO){
                if (me==master) cout << "  - Schreibe das Multigrid im Geomview-Format in ein File ... ";
                PrintGEO(pmg);
                if (me==master) cout << " OK\n";
            }

            CheckParMultiGrid(pmg,REF,master);

            DoMigration(pmg, LoadBal,C.coarseStrategy);

            movedCoarseNodes += LoadBal.GetMovedMultiNodes();

            if (C.printSize){
                if (me==master)
                    cout << "  - Verteilung der Elemente:\n";
                mg.SizeInfo(cout);
            }
            if (C.printPMG)
            {
                if (me==master) cout << "  - Schreibe Debug-Informationen in ein File ... ";
                PrintMG(pmg, MIG);
                if (me==master) cout << " OK\n";
            }

            CheckParMultiGrid(pmg,MIG,master);

            if (me==master && ref!=C.coarsedrop+C.coarseall-1) cout << line << endl;
        }

        if (me==master)
            cout << dline<< endl;

        if (C.printTime)
            Times.Print(cout);
        if (me==master)
            cout << "Moved Multinodes for refinement: " << movedRefNodes << endl
                    << "Moded Multinodes for coarsening: " << movedCoarseNodes << endl
                    << dline << endl << "Shuting down ...\n";

    }
    catch (DROPS::DROPSErrCL err) { cout << "In Assert gelaufen!\n";err.handle(); }
    return 0;
}


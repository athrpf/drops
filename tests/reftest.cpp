#include "geom/multigrid.h"
#include "geom/builder.h"
#include "out/output.h"
#include <fstream>
#include <cmath>
#include <cstdio>

void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    long ct= 0;
    
    for (DROPS::MultiGridCL::TetraIterator It(mg.GetAllTetraBegin(maxLevel)),
             ItEnd(mg.GetAllTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( It->IsUnrefined() && (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
        {
            It->SetRegRefMark(); ++ct;
        }
    }
    std::cerr << ct <<" Tetras wurden markiert\n";
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TetraIterator It(mg.GetAllTetraBegin(maxLevel)),
             ItEnd(mg.GetAllTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( It->IsUnrefined() && (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}


int main()
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.0;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    DROPS::MultiGridCL mg(brick);

//        std::cout << DROPS::DumpMGCL(mg) << std::endl << "hallo" << std::endl;
    mg.MakeConsistentNumbering();
    
//        std::cout << DROPS::DumpMGCL(mg);
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    mg.SizeInfo(std::cerr);
    for (DROPS::Uint i=0; i<2; ++i)
    {
        DROPS::MarkAll(mg);
        std::cerr << i+1 << ". Totalverfeinerung ------------------"<<std::endl;
        mg.Refine();
        std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
        mg.SizeInfo(std::cerr);
    }
//    std::cerr << DROPS::SanityMGOutCL(mg);
    for (DROPS::Uint i=0; i<4; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        DROPS::MarkAll(mg);
//        MarkSome(mg);
//        std::cerr << DROPS::SanityMGOutCL(mg);
//        std::cerr << DROPS::DumpMGCL(mg);
        std::cerr << i+1 << ". Tropfenverfeinerung ------------------"<<std::endl;
        mg.Refine();
        mg.SizeInfo(std::cerr);
        std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
//        DebugIt(&mg);
        char str[20];
        std::ofstream ofs;
        sprintf(str, "drop%i.off", i+1);
        ofs.open(str);
        ofs << DROPS::GeomMGOutCL(mg, -1, false, 0.5) << std::endl;
        ofs.close();
        sprintf(str, "bary%i.mg", i+1);
        ofs.open(str);
        for (DROPS::MultiGridCL::TriangTetraIteratorCL it= mg.GetTriangTetraBegin(), end= mg.GetTriangTetraEnd(); it!=end; ++it)
            ofs << GetBaryCenter(*it) << '\n';
        ofs.close();
        sprintf(str, "verts%i.mg", i+1);
        ofs.open(str);
        for (DROPS::MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(), end= mg.GetTriangVertexEnd(); it!=end; ++it)
            ofs << it->GetCoord() << '\n';
        ofs.close();
    }

//    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
//    mg.SizeInfo(std::cerr);

    int wait;
    std::cerr << "Press a key: " << std::flush;
    std::cin>>wait;

    for (DROPS::Uint i=0; i<6; ++i)
    {
        std::cerr << i+1 << ". Tropfenentfeinerung ------------------"<< std::endl;
        UnMarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        DROPS::UnMarkAll(mg);
        mg.Refine();
        std::cerr << DROPS::SanityMGOutCL(mg)  << std::endl;
        mg.SizeInfo(std::cerr);
    }
//    DROPS::UnMarkAll(mg);
//    for (DROPS::MultiGridCL::TetraIterator it= mg.GetAllTetraBegin(); it!=mg.GetAllTetraEnd(); ++it)
//        if (it->IsMarkedForRemovement()) std::cout << "MFR ";
//    std::cout << DROPS::DumpMGCL(mg);

    for (DROPS::Uint i=0; i<3; ++i)
    {
        DROPS::UnMarkAll(mg);
        std::cerr << i+1 << ". Totalentfeinerung ------------------"<<std::endl;
        mg.Refine();
        std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
        mg.SizeInfo(std::cerr);
    }

    {std::ofstream os("ttt2.off");
    os << DROPS::GeomMGOutCL(mg, -1, false) << std::endl;
    os.close();}
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

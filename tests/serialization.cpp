#include "misc/container.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include <fstream>

using namespace DROPS;
Uint rule = 0;

int main (int argc, char** argv)
{

    if (argc>1)
        rule = atoi(argv[1]);

    // Multigrid aufbauen
//    std::ifstream meshfile( "mesh.txt");
//    if (!meshfile) {
//        std::cout << "error while opening mesh file " << "mesh.txt" << "\n";
//        return 1;
//    }
//    DROPS::ReadMeshBuilderCL builder( meshfile);

    TetraBuilderCL builder(rule);
    MultiGridCL mg( builder);
//    MarkAll( mg);
//    mg.Refine();
//    MarkAll( mg);
//    mg.Refine();
    mg.SizeInfo( std::cout);
    std::cout << SanityMGOutCL(mg) << std::endl;

    MGSerializationCL serial (mg, "out-");
    serial.WriteMG();

//    std::ifstream meshfile2( "mesh.txt");
//    DROPS::ReadMeshBuilderCL builder2( meshfile2);
    TetraBuilderCL builder2(rule);
    FileBuilderCL buildmg("out-", &builder2);
    MultiGridCL mg2(buildmg);
    std::cout << "\n \n MultiGrid aus Datei gelesen\n \n";

    mg2.SizeInfo( std::cout);
    std::cout << SanityMGOutCL(mg2) << std::endl;

    MGSerializationCL serial2 (mg2, "neu-");
    serial2.WriteMG();
    return 0;
}

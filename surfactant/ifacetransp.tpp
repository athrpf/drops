/// \file
/// \brief Discretization for PDEs on an interface.

#include "levelset/levelset.h"
#include <cstring>

namespace DROPS {

template <class DiscVelSolT>
void SetupConvectionP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const DiscVelSolT& u)
{
    const IdxT num_unks= mat->RowIdx->NumUnknowns();
    MatrixBuilderCL m( &mat->Data, num_unks, num_unks);
    const Uint lvl= mat->GetRowLevel();
    IdxT num[4];

    std::cerr << "entering SetupConvectionP1: " << num_unks << " dof. ";

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> qp1[4];

    LocalP2CL<Point3DCL> u_loc;
    Point3DCL grad[4];

    double coup[4][4];
    double dummy;

    InterfacePatchCL patch;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *mat->RowIdx);
            P1DiscCL::GetGradients( grad, dummy, *it);
            u_loc.assign( *it, u);
            std::memset( coup, 0, 4*4*sizeof( double));

            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri)
                    SetupConvectionP1OnTriangle( &patch.GetBary( tri), patch.GetFuncDet( tri),
                        p1, qp1, u_loc, grad, coup);
            }

            for(int i= 0; i < 4; ++i) {// assemble row Numb[i]
                if (num[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (num[j] == NoIdx) continue;
                    m( num[i], num[j])+= coup[i][j]; // Order of indices is correct as the assemply of coup is adapted.
                }
            }
        }
    }
    m.Build();
    std::cerr << mat->Data.num_nonzeros() << " nonzeros in interface convection matrix!" << std::endl;
}


template <class DiscVelSolT>
void SetupMassDivP1 (const MultiGridCL& mg, MatDescCL* mat, const VecDescCL& ls, const DiscVelSolT& u)
{
    const IdxT num_unks= mat->RowIdx->NumUnknowns();
    MatrixBuilderCL m( &mat->Data, num_unks, num_unks);
    const Uint lvl= mat->GetRowLevel();
    IdxT num[4];

    std::cerr << "entering SetupMassDivP1: " << num_unks << " dof. ";

    LocalP1CL<> p1[4];
    p1[0][0]= p1[1][1]= p1[2][2]= p1[3][3]= 1.; // P1-Basis-Functions
    Quad5_2DCL<double> qp1[4];

    LocalP2CL<Point3DCL> u_loc;
    SMatrixCL<3,3> T;
    LocalP1CL<Point3DCL> gradrefp2[10], gradp2[10];
    P2DiscCL::GetGradientsOnRef( gradrefp2);

    double coup[4][4];
    double dummy;

    InterfacePatchCL patch;

    DROPS_FOR_TRIANG_CONST_TETRA( mg, lvl, it) {
        patch.Init( *it, ls);
        if (patch.Intersects()) { // We are at the phase boundary.
            GetLocalNumbP1NoBnd( num, *it, *mat->RowIdx);
            GetTrafoTr( T, dummy, *it);
            P2DiscCL::GetGradients( gradp2, gradrefp2, T);
            u_loc.assign( *it, u);
            std::memset( coup, 0, 4*4*sizeof( double));

            for (int ch= 0; ch < 8; ++ch) {
                patch.ComputeForChild( ch);
                for (int tri= 0; tri < patch.GetNumTriangles(); ++tri)
                    SetupMassDivP1OnTriangle( &patch.GetBary( tri), patch.GetFuncDet( tri),
                        p1, qp1, u_loc, gradp2, patch.GetNormal(), coup);
            }

            for(int i= 0; i < 4; ++i) {// assemble row Numb[i]
                if (num[i] == NoIdx) continue;
                for(int j= 0; j < 4; ++j) {
                    if (num[j] == NoIdx) continue;
                    m( num[i], num[j])+= coup[i][j];
                }
            }
        }
    }
    m.Build();
    std::cerr << mat->Data.num_nonzeros() << " nonzeros in mass-divergence matrix!" << std::endl;
}

} // end of namespace DROPS

//**************************************************************************
// File:    levelset.tpp                                                   *
// Content: levelset equation for two phase flow problems                  *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#include "num/discretize.h"
#include "levelset/fastmarch.h"
#include <fstream>

namespace DROPS
{

template<class DiscVelSolT>
void LevelsetP2CL::SetupSystem( const DiscVelSolT& vel)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.GetLevel(),
               idx= Phi.RowIdx->GetIdx();

    SparseMatBuilderCL<double> E(&_E, num_unks, num_unks), 
                               H(&_H, num_unks, num_unks);
    IdxT Numb[10];
    
    std::cerr << "entering SetupSystem: " << num_unks << " levelset unknowns. ";

    // fill value part of matrices
    Quad2CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad2CL<double> u_Grad[10]; // fuer u grad v_i
    SMatrixCL<3,3> T;
    double det, absdet, h_T;

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(_MG).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(_MG).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= fabs( det);
        h_T= pow( absdet, 1./3.);
        
        // collect some information about the edges and verts of the tetra
        // and save it in Numb and u_loc
        for(int i=0; i<4; ++i)
        {
            Numb[i]= sit->GetVertex(i)->Unknowns(idx);
	    u_loc[i]= vel.val( *sit->GetVertex(i));
        }
        for(int i=0; i<6; ++i)
        {
            Numb[i+4]= sit->GetEdge(i)->Unknowns(idx);
        }
        u_loc[4]= vel.val( *sit, 0.25, 0.25, 0.25);

        for(int i=0; i<10; ++i)
            u_Grad[i]= dot( u_loc, Grad[i]);
        
        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            for(int j=0; j<10; ++j)
            {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
                E( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j) * absdet
                                     + u_Grad[i].quadP2(j, absdet)*_SD*h_T; 
                
                // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                H( Numb[i], Numb[j])+= u_Grad[j].quadP2(i, absdet)
                                     + Quad2CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * _SD*h_T;
            }
    }
    
    E.Build();
    H.Build();
    std::cerr << _E.num_nonzeros() << " nonzeros in E, "
              << _H.num_nonzeros() << " nonzeros in H! " << std::endl;
}

} // end of namespace DROPS


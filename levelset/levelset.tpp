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
    const Uint lvl= Phi.GetLevel();

    SparseMatBuilderCL<double> E(&E_, num_unks, num_unks), 
                               H(&H_, num_unks, num_unks);
    IdxT Numb[10];
    
    std::cerr << "entering SetupSystem: " << num_unks << " levelset unknowns. ";

    // fill value part of matrices
    Quad2CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad2CL<double> u_Grad[10]; // fuer u grad v_i
    SMatrixCL<3,3> T;
    double det, absdet, h_T;

    P2DiscCL::GetGradientsOnRef( GradRef);
    
    for (MultiGridCL::const_TriangTetraIteratorCL sit=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(lvl), send=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd(lvl);
         sit!=send; ++sit)
    {
        GetTrafoTr( T, det, *sit);
        P2DiscCL::GetGradients( Grad, GradRef, T);
        absdet= std::fabs( det);
        h_T= std::pow( absdet, 1./3.);

        // save information about the edges and verts of the tetra in Numb
        GetLocalNumbP2NoBnd( Numb, *sit, *Phi.RowIdx);

        // save velocities inside tetra for quadrature in u_loc
        for(int i=0; i<4; ++i)
	    u_loc[i]= vel.val( *sit->GetVertex(i));
        u_loc[4]= vel.val( *sit, 0.25, 0.25, 0.25);

        for(int i=0; i<10; ++i)
            u_Grad[i]= dot( u_loc, Grad[i]);
        
        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            for(int j=0; j<10; ++j)
            {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
                E( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j) * absdet
                                     + u_Grad[i].quadP2(j, absdet)*SD_*h_T; 
                
                // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                H( Numb[i], Numb[j])+= u_Grad[j].quadP2(i, absdet)
                                     + Quad2CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * SD_*h_T;
            }
    }
    
    E.Build();
    H.Build();
    std::cerr << E_.num_nonzeros() << " nonzeros in E, "
              << H_.num_nonzeros() << " nonzeros in H! " << std::endl;
}

template<class ValueT>
ValueT InterfacePatchCL::quad( const LocalP2CL<ValueT>& f, double absdet, bool part /*, bool debug*/)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch_]);
    typedef BaryCoordCL* BaryPtrT;
    BaryPtrT BaryPtr[4];
    if (intersec_<3)
    { // cuts = Kind + leere Menge
//if (debug) std::cerr <<"Fall <3:\tKind + leere Menge\n";
        if ( part == (num_sign_[2]>0) )
        { // integriere ueber Kind
            for (int i=0; i<4; ++i)
                BaryPtr[i]= &BaryDoF_[data.Vertices[i]];
            return P1DiscCL::Quad( f, BaryPtr)*(absdet/8);
        }
        else
            return ValueT();
    }
    else if (intersec_==3)
    { // cuts = Tetra + Tetra-Stumpf
//if (debug) std::cerr <<"Fall 3:\tTetra + Tetra-Stumpf\n";
        int vertA= -1;  // cut-Tetra = APQR
        const int signA= num_sign_[0]==1 ? -1 : 1; // sign of vert A occurs only once
        for (int i=0; vertA==-1 && i<4; ++i)
            if (sign_[data.Vertices[i]]==signA) vertA= i;
        for (int i=0; i<3; ++i)
            BaryPtr[i]= &Bary_[i];
        BaryPtr[3]= &BaryDoF_[data.Vertices[vertA]];
        
        const double volFrac= VolFrac( BaryPtr);
        const ValueT quadTetra= P1DiscCL::Quad( f, BaryPtr)*(absdet*volFrac);
//if (debug) std::cerr << "vertA = " << vertA << "\tvolFrac = " << volFrac << "\t= 1 / " << 1/volFrac << std::endl;
        if ( part == (signA==1) )
            return quadTetra;
        else // Gesamt-Tetra - cut-Tetra
        {
            for (int i=0; i<4; ++i)
                BaryPtr[i]= &BaryDoF_[data.Vertices[i]];
            return P1DiscCL::Quad( f, BaryPtr)*(absdet/8) - quadTetra;
        }
    }
    else // intersec_==4
    { // cuts = 2x 5-Flaechner mit 6 Knoten. connectivity:     R---------S
      //                                                      /|        /|
      //                                                     A-|-------B |
      //                                                      \|        \|
      //                                                       P---------Q
//if (debug) std::cerr <<"Fall 4: 2x 5-Flaechner\t";
        int vertAB[2], // cut mit VZ==part = ABPQRS
            signAB= part ? 1 : -1;
        for (int i=0, k=0; i<4 && k<2; ++i)
            if (sign_[data.Vertices[i]]==signAB) vertAB[k++]= i;
        // connectivity AP automatisch erfuellt, check for connectivity AR
        const bool AR= vertAB[0]==VertOfEdge(Edge_[2],0) || vertAB[0]==VertOfEdge(Edge_[2],1);
//if (debug) if (!AR) std::cerr << "\nAR not connected!\n";

//if (debug) std::cerr << "vertA = " << vertAB[0] << "\tvertB = " << vertAB[1] << std::endl;
//if (debug) std::cerr << "PQRS on edges\t"; for (int i=0; i<4; ++i) std::cerr << Edge_[i] << "\t"; std::cerr << std::endl;
        // Integriere ueber Tetras ABPR, QBPR, QBSR    (bzw. mit vertauschten Rollen von Q/R)
        // ABPR    (bzw. ABPQ)
        BaryPtr[0]= &BaryDoF_[data.Vertices[vertAB[0]]];    
        BaryPtr[1]= &BaryDoF_[data.Vertices[vertAB[1]]]; 
        BaryPtr[2]= &Bary_[0];    BaryPtr[3]= &Bary_[AR ? 2 : 1]; 
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cerr << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        ValueT integral= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        // QBPR    (bzw. RBPQ)
        BaryPtr[0]= &Bary_[AR ? 1 : 2];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cerr << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        integral+= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        // QBSR    (bzw. RBSQ)
        BaryPtr[2]= &Bary_[3];
//if (debug) { const double volFrac= VolFrac(BaryPtr); std::cerr << "volFrac = " << volFrac << " = 1 / " << 1/volFrac << std::endl; }
        integral+= P1DiscCL::Quad( f, BaryPtr)*VolFrac(BaryPtr);
        
        return absdet*integral;
    }    
}

} // end of namespace DROPS


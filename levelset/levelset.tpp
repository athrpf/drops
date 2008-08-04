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
void LevelsetP2CL::GetInfo( double& maxGradPhi, double& Volume, Point3DCL& bary, Point3DCL& vel, const DiscVelSolT& velsol, Point3DCL& minCoord, Point3DCL& maxCoord) const
/** 
 * - \p maxGradPhi is the maximal 2-norm of the gradient of the level set function. This can be used as an indicator to decide 
 *   whether a reparametrization should be applied.
 * - \p Volume is the volume inside the approximate interface consisting of planar segments.
 * - \p bary is the barycenter of the droplet.
 * - \p vel is the velocity of the barycenter of the droplet.
 * - The entries of \p minCoord store the minimal x, y and z coordinates of the approximative interface, respectively. 
 * - The entries of \p maxCoord store the maximal x, y and z coordinates of the approximative interface, respectively.
 */ 
{
    Quad2CL<Point3DCL> Grad[10], GradRef[10];
    SMatrixCL<3,3> T;
    double det, absdet;
    InterfacePatchCL patch;

    P2DiscCL::GetGradientsOnRef( GradRef);
    maxGradPhi= -1.;
    Volume= 0.;
    bary[0]= bary[1]= bary[2]= 0;
    vel[0]= vel[1]= vel[2]= 0;
    minCoord[0]= minCoord[1]= minCoord[2]= 1e99;
    maxCoord[0]= maxCoord[1]= maxCoord[2]= -1e99;
    LocalP2CL<double> ones( 1.); 
    LocalP2CL<Point3DCL> Coord, Vel;
    
    for (MultiGridCL::const_TriangTetraIteratorCL it=const_cast<const MultiGridCL&>(MG_).GetTriangTetraBegin(), end=const_cast<const MultiGridCL&>(MG_).GetTriangTetraEnd();
        it!=end; ++it)
    {
        GetTrafoTr( T, det, *it);
        absdet= std::abs( det);
        P2DiscCL::GetGradients( Grad, GradRef, T); // Gradienten auf aktuellem Tetraeder
    
        patch.Init( *it, Phi);
        
        // compute maximal norm of grad Phi
        Quad2CL<Point3DCL> gradPhi;
        for (int v=0; v<10; ++v) // init gradPhi, Coord
        {
            gradPhi+= patch.GetPhi(v)*Grad[v];
            Coord[v]= v<4 ? it->GetVertex(v)->GetCoord() : GetBaryCenter( *it->GetEdge(v-4));
        }
        Vel.assign( *it, velsol);
        VectorCL normGrad( 5); 
        for (int v=0; v<5; ++v) // init normGrad
            normGrad[v]= norm( gradPhi[v]);
        const double maxNorm= normGrad.max();
        if (maxNorm > maxGradPhi) maxGradPhi= maxNorm;
    
        for (int ch=0; ch<8; ++ch)
        {
            // compute volume, barycenter and velocity
            patch.ComputeCutForChild(ch);
            Volume+= patch.quad( ones, absdet, false);
            bary+= patch.quad( Coord, absdet, false);
            vel+= patch.quad( Vel, absdet, false);

            // find minimal/maximal coordinates of interface
            if (!patch.ComputeForChild(ch)) // no patch for this child
                continue;
            for (Uint i=0; i<patch.GetNumPoints(); ++i)
            {
                const Point3DCL p= patch.GetPoint(i);
                for (int j=0; j<3; ++j)
                {
                    if (p[j] < minCoord[j]) minCoord[j]= p[j];
                    if (p[j] > maxCoord[j]) maxCoord[j]= p[j];
                }
            }
        }
    }
    bary/= Volume;    vel/= Volume;
} 


template<class DiscVelSolT>
void LevelsetP2CL::SetupSystem( const DiscVelSolT& vel)
// Sets up the stiffness matrices:
// E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
// H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
// where v_i, v_j denote the ansatz functions.
{
    const IdxT num_unks= Phi.RowIdx->NumUnknowns;
    const Uint lvl= Phi.GetLevel();

    SparseMatBuilderCL<double> bE(&E, num_unks, num_unks),
                               bH(&H, num_unks, num_unks);
    IdxT Numb[10];

    std::cerr << "entering SetupSystem: " << num_unks << " levelset unknowns. ";

    // fill value part of matrices
    Quad5CL<Point3DCL> Grad[10], GradRef[10], u_loc;
    Quad5CL<double> u_Grad[10]; // fuer u grad v_i
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
        u_loc.assign( *sit, vel);

        for(int i=0; i<10; ++i)
            u_Grad[i]= dot( u_loc, Grad[i]);

//        double maxV = 0.; // scaling of SD parameter (cf. master thesis of Rodolphe Prignitz)
//        const double maxV_limit= 1e-5;
//        for(int i=0; i<Quad5CL<>::NumNodesC; ++i)
//            maxV = std::max( maxV, u_loc[i].norm());
//        if( maxV < maxV_limit) maxV= maxV_limit; // no scaling for extremely small velocities
//        /// \todo fixed limit for maxV (maxV_limit), any better idea?
        double maxV= 1; // no scaling
        for(int i=0; i<10; ++i)    // assemble row Numb[i]
            for(int j=0; j<10; ++j)
            {
                // E is of mass matrix type:    E_ij = ( v_j       , v_i + SD * u grad v_i )
                bE( Numb[i], Numb[j])+= P2DiscCL::GetMass(i,j) * absdet
                                      + u_Grad[i].quadP2(j, absdet)*SD_/maxV*h_T;

                // H describes the convection:  H_ij = ( u grad v_j, v_i + SD * u grad v_i )
                bH( Numb[i], Numb[j])+= u_Grad[j].quadP2(i, absdet)
                                      + Quad5CL<>(u_Grad[i]*u_Grad[j]).quad( absdet) * SD_/maxV*h_T;
            }
    }

    bE.Build();
    bH.Build();
    std::cerr << E.num_nonzeros() << " nonzeros in E, "
              << H.num_nonzeros() << " nonzeros in H! " << std::endl;
}

inline double
InterfacePatchCL::EdgeIntersection (Uint v0, Uint v1, LocalP2CL<> & philoc)
{
    if (LinearEdgeIntersection) return philoc[v0]/(philoc[v0]-philoc[v1]);
    else {
        const double l0= philoc[v0], l1= philoc[v1],
            lm= philoc( AllEdgeBaryCenter_[v0][v1]); // Value of Phi in the edge-barycenter.
        // If l0*l1<0 the quadratic equation with p(0)=l0, p(1)=l1, p(1/2)=lm has exactly one root in (0,1).
        const double quadcoeff= 2.*l0 + 2.*l1 - 4.*lm, lincoeff= 4.*lm - 3.*l0 - l1;
        if ( std::fabs( quadcoeff) < std::fabs( lincoeff)*8.*std::numeric_limits<double>::epsilon()) // linear LS-function
            return l0/(l0 - l1);
        const double rt= std::sqrt( std::pow(4.*lm - (l0 + l1), 2) - 4.*l0*l1);
        const double x0= (-lincoeff - rt)/(2.*quadcoeff),
                     x1= (-lincoeff + rt)/(2.*quadcoeff);
        Assert( (0 < x0 && x0 < 1.) || (0 < x1 && x1 < 1.),
            "InterfacePatchCL::EdgeIntersection: Excessive roundoff-error with quadratic level-set-function",
            DebugNumericC);
        return (0. < x0 && x0 < 1.) ? x0 : x1;
    }
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

template<class ValueT>
void InterfacePatchCL::quadBothParts( ValueT& int_pos, ValueT& int_neg, const LocalP2CL<ValueT>& f, double absdet)
{
    const ChildDataCL data= GetChildData( RegRef_.Children[ch_]);
    typedef BaryCoordCL* BaryPtrT;
    BaryPtrT BaryPtr[4], chTetra[4];

    for (int i=0; i<4; ++i)
        chTetra[i]= &BaryDoF_[data.Vertices[i]];
    const ValueT quadChild= P1DiscCL::Quad( f, chTetra)*(absdet/8);
    
    if (intersec_<3)
    { // cuts = Kind + leere Menge
//if (debug) std::cerr <<"Fall <3:\tKind + leere Menge\n";
        if ( num_sign_[2]>0 )
        {
            int_pos= quadChild;    int_neg= ValueT();
        }
        else
        {
            int_neg= quadChild;    int_pos= ValueT();
        }
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
        if ( signA==1 )
        {
            int_pos= quadTetra;   int_neg= quadChild - quadTetra;
        }
        else
        {
            int_neg= quadTetra;   int_pos= quadChild - quadTetra;
        }
    }
    else // intersec_==4
    { // cuts = 2x 5-Flaechner mit 6 Knoten. connectivity:     R---------S
      //                                                      /|        /|
      //                                                     A-|-------B |
      //                                                      \|        \|
      //                                                       P---------Q
//if (debug) std::cerr <<"Fall 4: 2x 5-Flaechner\t";
        int vertAB[2], // cut mit VZ == + = ABPQRS
            signAB= 1;
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

        int_pos= absdet*integral;    int_neg= quadChild - int_pos;
    }
}

template <class It>
  std::pair<double, double>
  h_interface (It begin, It end, const VecDescCL& ls)
{
    double hmin= 1e99, hmax= -1.0; //, hmean= 0.;
    // size_t num= 0;
    // EdgeCL* emax= 0;

    for (; begin != end; ++begin)
        if (   (ls.Data[begin->GetVertex( 0)->Unknowns( ls.RowIdx->GetIdx())]
                *ls.Data[begin->Unknowns( ls.RowIdx->GetIdx())] <= 0.)
            || (ls.Data[begin->GetVertex( 1)->Unknowns( ls.RowIdx->GetIdx())]
                *ls.Data[begin->Unknowns( ls.RowIdx->GetIdx())] <= 0.)) {
            const double h= (begin->GetVertex( 0)->GetCoord() - begin->GetVertex( 1)->GetCoord()).norm();
            hmin= std::min( hmin, h);
            // if (h > hmax) emax= &*begin;
            hmax= std::max( hmax, h);
            // hmean+= h;
            // ++num;
        }
    // std::cerr << "mean : " << hmean/num << '\n';
    // emax->DebugInfo( std::cerr);
    // emax->GetVertex( 0)->DebugInfo( std::cerr);
    // emax->GetVertex( 1)->DebugInfo( std::cerr);

    return std::make_pair( hmin, hmax);
}

} // end of namespace DROPS



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

namespace DROPS
{

void SetupPiPj(LocalP2CL<>pipj[4][4])
{    
    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            pipj[j][i][EdgeByVert( i, j) + 4]= 0.25;
            pipj[i][j][EdgeByVert( j, i) + 4]= 0.25;
        }
        pipj[i][i][i]= 1.;
        for (int vert= 0; vert < 3; ++vert)
            pipj[i][i][EdgeByVert( i, VertOfFace( i, vert)) + 4]= 0.25;
    }
}

bool IsRegBaryCoord(const SArrayCL<BaryCoordCL,4>& T)
{
    for(int i= 0; i < 4; ++i)  
        for(int j= 0; j < 4; ++j) {
            if (isnan(T[i][j])|| isinf(T[i][j])|| T[i][j] >1.|| T[i][j] <0.) {
                std::cout << "Irregular coordinate!\n";
                return false;
        }
    }
    return true;
}
//====================================================
//
// Setup routines in each tetra
//
//====================================================
void SetupLocalRhs( double locf[4], const Quad5CL<> rhs,  const Quad5CL<> p[4], const double absdet) 
{
    for(int i= 0; i < 4; ++i)
        locf[i]= Quad5CL<>(rhs*p[i]).quad( absdet);
}

void SetupLocalTwoPhaseRhs(const TetraCL& t, InterfaceTetraCL& cut, const LocalP1CL<> p1[4],  double f_n[4], double f_p[4], 
    const instat_scalar_fun_ptr f, const double absdet, const double time) 
{
    cut.ComputeSubTets();
    for(int i= 0; i < 4; ++i) 
       f_n[i] = f_p[i] = 0.;
           
    for (Uint k=0; k< cut.GetNumTetra(); ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        if (!IsRegBaryCoord(T))  continue;
        double Vol = absdet*VolFrac(T);
        if (isnan(Vol) || isinf(Vol)){
            std::cout<<" SetupLocalTwoPhaseRhs: M Support of XFEM is too small.\t";
            continue;
        }
        BaryCoordCL* nodes;
        nodes = Quad5CL<>::TransformNodes(T);
        Quad5CL<> qrhs(t, f, time, nodes);  
        
        for(Uint i= 0; i < 4; ++i) {
            Quad5CL<> qp1(p1[i], nodes);
            double tmp = Quad5CL<>(qrhs*qp1).quad(Vol);
            if (isnan(tmp) || isinf(tmp)) {
                for( Uint j=0; j<4; ++j) {
                    f_p[j]= 0.;
                    f_n[j]= 0.;
                }    
                break;
            }
            if (k>= cut.GetNumNegTetra())
                f_p[i]+= tmp;
            else
                f_n[i]+= tmp;
        }
        if (nodes) delete nodes;
        nodes=0;
    }     
}

// Couplings between the XFEM basis functions (at the same time step or different time steps) when the tetra is cut by only one interface.
void SetupLocalOneInterfaceSystem( InterfaceTetraCL& cut, double M_n[4][4], double M_p[4][4], double A_n[4][4], double A_p[4][4], 
    double C_n[4][4], double C_p[4][4], const double absdet, const double D[2], const double H, const LocalP2CL<Point3DCL> lp2_u, 
    const LocalP2CL<> pipj[4][4], const LocalP1CL<> p1[4], const SMatrixCL<3,4> G, const SMatrixCL<4,4> GTG, bool sign[4], 
    bool jumps, bool pPart)
{
    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumTetra(); //# of subtetras
    // D, H are constant if tetra is cut by new interface
    double d = pPart ? D[0] : D[1]/H,
           h = pPart ? 1. : 1./H;
    
    std::memset( M_n, 0, 4*4*sizeof(double));
    std::memset( A_n, 0, 4*4*sizeof(double));
    std::memset( C_n, 0, 4*4*sizeof(double));
    std::memset( M_p, 0, 4*4*sizeof(double));
    std::memset( A_p, 0, 4*4*sizeof(double));
    std::memset( C_p, 0, 4*4*sizeof(double));

    for(int i= 0; i < 4; ++i) 
        sign[i]= (cut.GetSign(i) == 1);
    
    for (Uint k=0; k< NumTets; ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        if (!IsRegBaryCoord(T))  continue;
        double Vol = absdet*VolFrac(T);
        if (isnan(Vol)|| isinf(Vol)){
            std::cout<<"Vol " <<VolFrac(T)<<"\n";
            std::cout<<"SetupLocalOneInterfaceSystem: Support of XFEM is too small.\t";
            continue;
        }
        BaryCoordCL* nodes;
        nodes = Quad3CL<>::TransformNodes(T);
        Quad3CL<Point3DCL> qu(lp2_u, nodes);
        bool irreg = false; 
        for(int i= 0; i < 4; ++i) {
            Quad3CL<> qp1(p1[i], nodes);
            for(int j= 0; j < 4; ++j) {
                Quad3CL<> qM(pipj[j][i],nodes),
                          qC(dot( qu, Quad3CL<Point3DCL>(G.col( j)))*qp1); // G is constant in t
                double iM = qM.quad(Vol), iC = qC.quad(Vol), iA = Vol* GTG( j, i)/6.; 
                if (isnan(iM)|| isinf(iM)||isnan(iA)|| isinf(iA)||isnan(iC)|| isinf(iC)) {
                    std::memset( M_n, 0, 4*4*sizeof(double));
                    std::memset( A_n, 0, 4*4*sizeof(double));
                    std::memset( C_n, 0, 4*4*sizeof(double));
                    std::memset( M_p, 0, 4*4*sizeof(double));
                    std::memset( A_p, 0, 4*4*sizeof(double));
                    std::memset( C_p, 0, 4*4*sizeof(double));
                    irreg = true;
                    break;
                }
       
                // D and H are piecewise constants in T. When
                // - Compute the integrals between XFEM basis functions at the same time step
                // - Compute the integrals between new FEM basis functions with old XFEM
                // basis functions (wrt old interface) -> (FEM function, old XFEM).
                if(jumps){
                    if (k<cut.GetNumNegTetra()){
                        M_n[j][i]+= iM/H;
                        A_n[j][i]+= iA*D[1]/H;
                        C_n[j][i]+= iC/H;
                    }
                    else{
                        M_p[j][i]+= iM ;
                        A_p[j][i]+= iA*D[0];
                        C_p[j][i]+= iC;
                    }
                }
                // D and H are constant.
                // Compute (new XFEM test function, FEM) at old time step. Tetra is cut only by new interface.
                else{  
                    if (k<cut.GetNumNegTetra()){
                        M_n[j][i]+= iM*h;
                        A_n[j][i]+= iA*d;
                        C_n[j][i]+= iC*h;
                    }
                    else{
                        M_p[j][i]+= iM*h;
                        A_p[j][i]+= iA*d;
                        C_p[j][i]+= iC*h;
                    }    
                }
            }
            if (irreg) break;
        }
        if(nodes) delete nodes;
        nodes=0;
    }
}

void SetupLocalOnePhaseMassMatrix( double locM[4][4], const double absdet, const double H, bool pPart)
{
    const double h = pPart ? 1. : 1./H;
    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            locM[j][i]= h*P1DiscCL::GetMass( i, j)*absdet;
            locM[i][j]= locM[j][i];
        }
        locM[i][i]= h*P1DiscCL::GetMass( i, i)*absdet;
    }
}



/// compute the mass matrix for the case that the tetrahedron is cut ONLY for the old OR new time. \n
/// computes the off-diagonal block of \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ M_{1,2}[i,j] = \int h \phi_j^{new} \phi_i^{old} \f$ \n
/// or
/// \f$ M_{2,1}[i,j] = \int h \phi_j^{old} \phi_i^{new} \f$ \n
void SetupLocalOneInterfaceMassMatrix( InterfaceTetraCL& cut, double M_n[4][4], double M_p[4][4], 
    const double absdet, const double H, const LocalP2CL<> pipj[4][4], bool sign[4], bool jumps, bool pPart)
{
    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumTetra(); /// # of subtetras
    // D, H are constant if T \cap Gamma_old is empty
    double h = pPart ? 1. : 1./H;
    std::memset( M_n,0, 4*4*sizeof(double));
    std::memset( M_p,0, 4*4*sizeof(double));
    
    for(int i= 0; i < 4; ++i) {
        sign[i]= (cut.GetSign(i) == 1);  
    }        

    for (Uint k=0; k< NumTets; ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        if (!IsRegBaryCoord(T)) continue;
        
        double Vol = absdet*VolFrac(T);
        if (isnan(Vol)|| isinf(Vol)){
            std::cout<<"Vol " <<VolFrac(T)<<"\n";
            std::cout<<" Support of XFEM is too small.\t";
            continue;
        }
        BaryCoordCL* nodes;
        nodes = Quad3CL<>::TransformNodes(T);
        bool irreg=false;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < 4; ++j) {
                Quad3CL<> qM(pipj[j][i],nodes);
                         
                double iM = qM.quad(Vol);
                if (isnan(iM)|| isinf(iM)) {
                    std::memset( M_n,0, 4*4*sizeof(double));
                    std::memset( M_p,0, 4*4*sizeof(double));
                    irreg=true;
                    break; //leave the j-loop 
                }

                // D and H are piecewise constants in T. When
                // - Compute the integrals between XFEM basis functions at the same time step
                // - Compute the integrals between new FEM basis functions with old XFEM
                // basis functions (wrt old interface) -> (FEM function, old XFEM).
                if(jumps){
                    if (k<cut.GetNumNegTetra())
                        M_n[j][i]+= iM/H;
                    else
                        M_p[j][i]+= iM ;
                }
                // D and H are constant.
                // Compute (new XFEM test function, FEM) at old time step. Tetra is cut by only new interface.
                else{  
                    if (k<cut.GetNumNegTetra())
                        M_n[j][i]+= iM*h;
                    else
                        M_p[j][i]+= iM*h;
                }
            }
            if (irreg) break;  //leave the i-loop
        }
        if (nodes) delete nodes;
        nodes=0;
    }
}

/// compute the following element matrices for an element which is completely in one phase only: \n
/// \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ A(u,v) = \int d \nabla u \nabla v \f$ with d the diffusion coefficient times the inverse of the henry coefficient of the domain \n
/// \f$ C(u,v) = \int (b \cdot \nabla u) v \f$ with b the local velocity times the inverse of the henry coefficient of the domain \n
void SetupLocalOnePhaseSystem( double locM[4][4], double locA[4][4], double locC[4][4], const Quad3CL<Point3DCL> q3_u, 
    const Quad3CL<> p[4], const SMatrixCL<3,4> G, const SMatrixCL<4,4> GTG, const double absdet, const double D[2], 
    const double H, bool pPart)
{
    const double d = pPart ? D[0] : D[1]/H;
    const double hw = pPart ? 1. : 1./H;
    for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < i; ++j) {
            locM[j][i]= hw*P1DiscCL::GetMass( i, j)*absdet;
            locM[i][j]= locM[j][i];
            locA[j][i]= d*GTG( i, j)/6.0*absdet;
            locA[i][j]= locA[j][i];
            locC[j][i]= hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( j)))*p[i]).quad( absdet);
            locC[i][j]= hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( i)))*p[j]).quad( absdet);
        }
        locM[i][i]= hw*P1DiscCL::GetMass( i, i)*absdet;
        locA[i][i]= d*GTG( i, i)/6.0*absdet;
        locC[i][i]= hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( i)))*p[i]).quad( absdet);
    }
}

/// compute the mass matrix for the case that the tetrahedron is cut for the old and new time. \n
/// computes \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ M_{2,1}[i,j] = \int h \phi_j^{old} \phi_i^{new} \f$ \n
/// \f$ M_{2,2}[i,j] = \int h \phi_j^{old} \phi_i^{old} \f$ \n
void SetupLocalTwoInterfacesMassMatrix( InterfaceTetraCL& cut, InterfaceTetraCL& oldcut, 
    double M22[4][4], double M21[4][4], const double absdet, const double H, LocalP2CL<>& lp2_oldlset,
    const LocalP2CL<> pipj[4][4])
{
    bool sign[4], oldsign[4];
    std::memset( M22, 0, 4*4*sizeof(double));
    std::memset( M21, 0, 4*4*sizeof(double));
    
    for(int i= 0; i < 4; ++i) {
        sign[i]= (cut.GetSign(i) == 1);
        oldsign[i]= (oldcut.GetSign(i) == 1);
    }        

    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumTetra(); /// # of subtetras

    for (Uint k=0; k< NumTets; ++k){
        bool Tk= (k>=cut.GetNumNegTetra());    // Tk in Omega_new_+?
        const SArrayCL<BaryCoordCL,4>& T =  cut.GetTetra(k);
        if (!IsRegBaryCoord(T)) continue;
        
        InterfaceTetraCL suboldcut;
        suboldcut.Init(T, lp2_oldlset,0.);
        double VolT = absdet*VolFrac(T);
        bool nocut= !suboldcut.Intersects();
        suboldcut.ComputeSubTets();
        Uint NumOldTets= suboldcut.GetNumTetra();
        if (nocut){  
            bool pPart = (suboldcut.GetSign( 0) == 1);
            double h = pPart ? 1. : 1./H;

            BaryCoordCL *nodes1;
            nodes1 = Quad3CL<>::TransformNodes(T); //nodes on physical domain
            bool irreg=false;
            for(int i= 0; i < 4; ++i) {
                if (Tk == sign[i]) continue; // supp(Phi_i^{Gamma_new}) \cap Tk is empty
                for(int j= 0; j < 4; ++j) {
                    Quad3CL<> qM(pipj[i][j],nodes1); 
                    double iM = qM.quad(VolT)*h;
                    if (isnan(iM)|| isinf(iM)) {
                    ///> TODO: If a local value in a child tetrahedron is irregular, ignore this tetra
                    ///> The contribution of other tetra must be preserved
                    ///> For example tmpM21[4][4]
//                         std::memset( M21,0, 4*4*sizeof(double));
//                         std::memset( M22,0, 4*4*sizeof(double));
                        irreg=true;
                        break; //leave the j-loop 
                    }
                    // int_(Phi_i_Gamma_new * Phi_j) in Tk
                    M21[j][i]+= sign[i]? -iM: iM;

                    if (pPart== oldsign[j]) continue; // supp(Phi_j^{Gamma_old}) \cap Tk is empty

                    M22[j][i]+= (sign[i]==oldsign[j]) ? iM : -iM ;
                }
                if (irreg) break;
            }
            if(nodes1) delete nodes1;
            nodes1=0;
            continue;
        }
           
        for (Uint m=0; m< NumOldTets; ++m){
            bool Tkm= (m>=suboldcut.GetNumNegTetra());  // Tkm in Omega_old_+?
            double h = Tkm ? 1. : 1./H;
            const SArrayCL<BaryCoordCL,4>& Tc =  suboldcut.GetTetra(m);
            if (!IsRegBaryCoord(Tc)) continue;
            double Vol = absdet*VolFrac(Tc);
            BaryCoordCL *nodes2;
            nodes2 = Quad3CL<>::TransformNodes(Tc);
            bool irreg=false;
            for(int i= 0; i < 4; ++i) {
                if (Tk == sign[i]) continue; // supp(Phi_i^{Gamma_new}) \cap Tk is empty
                for(int j= 0; j < 4; ++j) {
                    Quad3CL<> qM(pipj[i][j],nodes2); // G is constant in t
                    double iM = qM.quad(Vol)*h;
                    if (isnan(iM)|| isinf(iM)) {
                        ///> TODO: If a local value in a child tetrahedron is irregular, ignore this tetra
                        ///> The contribution of other tetra must be preserved
                        ///> For example tmpM21[4][4]
//                         std::memset( M21,0, 4*4*sizeof(double));
//                         std::memset( M22,0, 4*4*sizeof(double));
                        irreg=true;
                        break; //leave the j-loop 
                    }
                    M21[j][i]+= sign[i]? -iM: iM;
                    if (Tkm == oldsign[j]) continue;// supp(Phi_j^{Gamma_old}) \cap Tkm is empty

                    M22[j][i]+= (sign[i]==oldsign[j]) ? iM : -iM ;
                }
                if (irreg) break;
            }
            if (nodes2) delete nodes2;
            nodes2=0;
       }   
    }
}

// void SetupLocalTwoInterfacesSystem( InterfaceTetraCL& cut, InterfaceTetraCL& oldcut, 
//     double M22[4][4], double A22[4][4], double C22[4][4], double M21[4][4], double A21[4][4], double C21[4][4],
//     const double absdet, const double D[2], const double H, const LocalP2CL<Point3DCL> lp2_u, LocalP2CL<>& lp2_oldlset,
//     const LocalP2CL<> pipj[4][4], const LocalP1CL<> p1[4], const SMatrixCL<3,4> G,  const SMatrixCL<4,4> GTG) 
// {
//     bool sign[4], oldsign[4];
//     std::memset( M22,0, 4*4*sizeof(double));
//     std::memset( A22,0, 4*4*sizeof(double));
//     std::memset( C22,0, 4*4*sizeof(double));
//     std::memset( M21,0, 4*4*sizeof(double));
//     std::memset( A21,0, 4*4*sizeof(double));
//     std::memset( C21,0, 4*4*sizeof(double));
// 
//     for(int i= 0; i < 4; ++i) {
//         sign[i]= (cut.GetSign(i) == 1); 
//         oldsign[i]= (oldcut.GetSign(i) == 1);
//     }        
//     
//     cut.ComputeSubTets();
//     Uint NumTets=cut.GetNumTetra(); // # of subtetras
// 
//     for (Uint k=0; k< NumTets; ++k){
//         bool Tk= (k>=cut.GetNumNegTetra());    // Tk in Omega_new_+?
//         const SArrayCL<BaryCoordCL,4>& T =  cut.GetTetra(k);
//         if (!IsRegBaryCoord(T)) continue;
//         
//         InterfaceTetraCL suboldcut;
//         suboldcut.Init(T, lp2_oldlset);
//         double VolT = absdet*VolFrac(T);
//         if (isnan(VolT)|| isinf(VolT)){
//             std::cout<<"SetupLocalTwoInterfaceSystem: T Support of XFEM is too small.\t";
//             continue;
//         }
//         
//         bool nocut= !suboldcut.Intersects();
//         suboldcut.ComputeSubTets();
//         Uint NumOldTets= suboldcut.GetNumTetra();//std::cout<< NumOldTets <<"\t";
//         if (nocut ){  
//             bool pPart = suboldcut.GetSign( 0) == 1;
//             double d = pPart ? D[0] : D[1]/H;
//             double h = pPart ? 1. : 1./H;
// 
//             BaryCoordCL *nodes1;
//             nodes1 = Quad3CL<>::TransformNodes(T);
//             Quad3CL<Point3DCL> qu(lp2_u, nodes1);
//             bool irreg=false;
//             for(int i= 0; i < 4; ++i){
//                 if (Tk == sign[i]) continue; // supp(Phi_i^{Gamma_new}) \cap Tk is empty
//                 Quad3CL<> qp1(p1[i], nodes1);
//                 for(int j= 0; j < 4; ++j){
//                     Quad3CL<> qM(pipj[i][j],nodes1),
//                               qC(dot( qu, Quad3CL<Point3DCL>(G.col( j)))*qp1); // G is constant in t
//                     double iM = qM.quad(VolT)*h, iC = qC.quad(VolT)*h, iA = VolT* GTG( j, i)/6.*d;
//                     if (isnan(iM)|| isinf(iM)||isnan(iA)|| isinf(iA)||isnan(iC)|| isinf(iC)) {
//                         ///> TODO: If a local value in a child tetrahedron is irregular, ignore this tetra
//                         ///> The contribution of other tetra must be preserved
//                         ///> For example tmpM22[4][4]
// //                         std::memset( M22,0, 4*4*sizeof(double));
// //                         std::memset( A22,0, 4*4*sizeof(double));
// //                         std::memset( C22,0, 4*4*sizeof(double));
// //                         std::memset( M21,0, 4*4*sizeof(double));
// //                         std::memset( A21,0, 4*4*sizeof(double));
// //                         std::memset( C21,0, 4*4*sizeof(double));
//                         irreg=true;
//                         break; //leave the j-loop 
//                     }
//                     // int_(Phi_i_Gamma_new * Phi_j) in Tk
//                     M21[j][i]+= sign[i]? -iM: iM;
//                     A21[j][i]+= sign[i]? -iA: iA;
//                     C21[j][i]+= sign[i]? -iC: iC;
// 
//                     if (pPart== oldsign[j]) continue; // supp(Phi_j^{Gamma_old}) \cap Tk is empty
// 
//                     M22[j][i]+= (sign[i]==oldsign[j]) ? iM : -iM ;
//                     A22[j][i]+= (sign[i]==oldsign[j]) ? iA : -iA;
//                     C22[j][i]+= (sign[i]==oldsign[j]) ? iC : -iC;
//                 }
//                 if (irreg) break;
//             }
//             if (nodes1) delete nodes1;
//             nodes1=0;
//             continue;
//         }
//             
//         for (Uint m=0; m< NumOldTets; ++m){
//             bool Tkm= (m>=suboldcut.GetNumNegTetra());  // Tkm in Omega_old_+?
//             double d= Tkm ? D[0] : D[1]/H;
//             double h = Tkm ? 1. : 1./H;
//             const SArrayCL<BaryCoordCL,4>& Tc =  suboldcut.GetTetra(m);
//             if (!IsRegBaryCoord(Tc)) continue;
//             double Vol = absdet*VolFrac(Tc);
//             if (isnan(Vol)|| isinf(Vol)){
//                 std::cout<<" SetupLocalTwoInterfaceSystem: Tc Support of XFEM is too small.\t";
//                 continue;
//             }
// 
//             BaryCoordCL *nodes2;
//             nodes2 = Quad3CL<>::TransformNodes(Tc);
//             Quad3CL<Point3DCL> qu(lp2_u, nodes2);
// 
//             bool irreg=false;
//             for(int i= 0; i < 4; ++i) {
//                 if (Tk == sign[i]) continue; // supp(Phi_i^{Gamma_new}) \cap Tk is empty
//                 Quad3CL<> qp1(p1[i], nodes2);
//                 for(int j= 0; j < 4; ++j) {
//                     Quad3CL<> qM(pipj[i][j],nodes2),
//                               qC(dot( qu, Quad3CL<Point3DCL>(G.col( j)))*qp1); // G is constant in t
//                     
//                     double iM = qM.quad(Vol)*h, iC = qC.quad(Vol)*h, iA = Vol* GTG( j, i)/6.*d;
//                     if (isnan(iM)|| isinf(iM)||isnan(iA)|| isinf(iA)||isnan(iC)|| isinf(iC)) {
//                         ///> TODO: If a local value in a child tetrahedron is irregular, ignore this tetra
//                         ///> The contribution of other tetra must be preserved
//                         ///> For example tmpM22[4][4]
// //                         std::memset( M22,0, 4*4*sizeof(double));
// //                         std::memset( A22,0, 4*4*sizeof(double));
// //                         std::memset( C22,0, 4*4*sizeof(double));
// //                         std::memset( M21,0, 4*4*sizeof(double));
// //                         std::memset( A21,0, 4*4*sizeof(double));
// //                         std::memset( C21,0, 4*4*sizeof(double));
//                         irreg=true;
//                         break; //leave the j-loop 
//                     }
//                     M21[j][i]+= sign[i]? -iM: iM;
//                     A21[j][i]+= sign[i]? -iA: iA;
//                     C21[j][i]+= sign[i]? -iC: iC;
//                     if (Tkm == oldsign[j]) continue;// supp(Phi_j^{Gamma_old}) \cap Tkm is empty
//                     M22[j][i]+= (sign[i]==oldsign[j]) ? iM : -iM ;
//                     A22[j][i]+= (sign[i]==oldsign[j]) ? iA : -iA;
//                     C22[j][i]+= (sign[i]==oldsign[j]) ? iC : -iC;
// 
//                 }
//                 if (irreg) break;
//             }
//             if(nodes2) delete nodes2;
//             nodes2=0;
//        }   
//     }
// }


///Sets up the Nitsche System, i.e. the interface integrals for the formulation where H*u is the unknown!
void SetupLocalNitscheSystem( const BaryCoordCL * const p, const ExtIdxDescCL& Xidx, 
    Quad5_2DCL<Point3DCL> n, Point3DCL G[4], LocalNumbP1CL ln, MatrixBuilderCL& A, const double det, const double D[2], 
    const double H, const double kappa[2], const double lambda, const double h, const int sign[4])
{
    if (isnan(det)|| isinf(det)){
        std::cout<<" Support of XFEM function is too small.\t";
        return;
    }
    Quad5_2DCL<>    p1[4], Grad1_n[4];   // P1-Hat-Functions and their gradients...
    P1DiscCL::GetP1Basis( p1, p);
    for (int v=0; v<4; ++v)
        Grad1_n[v]=dot(G[v], n);

    double D_avr= D[0]*kappa[0]+ D[1]*kappa[1]/H, d[4];
    for(int i= 0; i < 4; ++i){
        d[i]= (sign[i]==1)? -D[1]*kappa[1]/H: D[0]*kappa[0];
    }
    for(int i= 0; i < 4; ++i)
        if(ln.WithUnknowns(i)){
            const IdxT xidx_i= Xidx[ln.num[i]];
            for(int j= 0; j < 4; ++j)
                if(ln.WithUnknowns(j)){
                const IdxT xidx_j= Xidx[ln.num[j]];
                if (xidx_j!=NoIdx ){
                    A( ln.num[i], xidx_j)-= Quad5_2DCL<>(D_avr*p1[j]*Grad1_n[i]).quad(det);
                }
                if (xidx_i!=NoIdx){
                    A( xidx_i, ln.num[j])-= Quad5_2DCL<>(D_avr*p1[i]*Grad1_n[j]).quad(det);
                }   
                if (xidx_i!=NoIdx && xidx_j!=NoIdx ){
                    A( xidx_i, xidx_j)-=  Quad5_2DCL<>(d[i]*p1[j]*Grad1_n[i] + d[j]*p1[i]*Grad1_n[j]
                    - lambda/h*p1[i]*p1[j]).quad(det);
                }
             }
        }
}

}//end of namespace DROPS

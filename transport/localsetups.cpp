
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

/// compute element vector for integrator: \n
/// \f$ f(v) = \int_T f v dx \f$ (for unstabilized fem)\n 
/// \f$ f(v) = \int_T f (v+\delta_T vel \cdot \nabla v) dx = \int_T f \tilde{v} dx \f$ (for sdfem)\n 
void ComputeRhsElementVector( double locf[4], LocalConvDiffReacCoefficients& local_cdcoef, TransformedP1FiniteElement& transformedfel)
{
    if (transformedfel.stabilized()){
      StabilizedTransformedP1FiniteElement& stabfe = static_cast<StabilizedTransformedP1FiniteElement&>(transformedfel);
      for(int i= 0; i < 4; ++i)
          locf[i]= Quad3CL<>(local_cdcoef.GetSourceAsQuad3()*stabfe.GetTestShapeAsQuad3CL(i)).quad(transformedfel.GetAbsDeterminant());
    }
    else
    {
      for(int i= 0; i < 4; ++i)
          locf[i]= Quad3CL<>(local_cdcoef.GetSourceAsQuad3()*transformedfel.GetGridfunctions().GetShapeAsQuad3(i)).quad(transformedfel.GetAbsDeterminant());
    }
}

void SetupLocalTwoPhaseRhs(TransformedP1FiniteElement& transformedfel, InterfaceTetraCL& cut, 
    ConvDiffElementVectors& elvecs, LocalConvDiffReacCoefficients& local_coefs) 
{
    cut.ComputeSubTets();
    
    elvecs.ResetSigned();//CL <--- kann weg

    StabilizedTransformedP1FiniteElement* stabfe = NULL;
    if (transformedfel.stabilized())
      stabfe = static_cast<StabilizedTransformedP1FiniteElement*>(&transformedfel);


    for (Uint k=0; k< cut.GetNumTetra(); ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        if (!IsRegBaryCoord(T))  continue;
        double Vol = transformedfel.GetAbsDeterminant()*VolFrac(T);
        if (isnan(Vol) || isinf(Vol)){
            std::cout<<" SetupLocalTwoPhaseRhs: M Support of XFEM is too small.\t";
            continue;
        }
        transformedfel.SetSubTetra(T);
        Quad3CL<> qrhs(transformedfel.GetTetra(), local_coefs.GetSourceAsFuntionPointer(), local_coefs.GetTime(), transformedfel.GetNodes());  
        
        for(Uint i= 0; i < 4; ++i) {
            double tmp;
            if (stabfe)
              tmp = Quad3CL<>(qrhs*stabfe->GetTestShapeAsQuad3CL(i)).quad(Vol);
            else
              tmp = Quad3CL<>(qrhs*transformedfel.GetBaseShapeAsQuad3CL(i)).quad(Vol);
            if (isnan(tmp) || isinf(tmp)) {
                for( Uint j=0; j<4; ++j) {
                    elvecs.f_p[j]= 0.;
                    elvecs.f_n[j]= 0.;
                }    
                break;
            }
            if (k>= cut.GetNumNegTetra())
                elvecs.f_p[i]+= tmp;
            else
                elvecs.f_n[i]+= tmp;
        }
    }     
}


/// compute the following element matrices for an element which is completely in one phase only: \n
/// \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ A(u,v) = \int d \nabla u \nabla v \f$ with d the diffusion coefficient times the inverse of the henry coefficient of the domain \n
/// \f$ C(u,v) = \int (b \cdot \nabla u) v \f$ with b the local velocity times the inverse of the henry coefficient of the domain \n
void SetupLocalOnePhaseSystem(TransformedP1FiniteElement& transformedfel, ConvDiffElementMatrices & elmats, LocalConvDiffReacCoefficients& local_coefs,	bool pPart)
{
    const SMatrixCL<3,4> G = transformedfel.GetDShape();
    const SMatrixCL<4,4> GTG = transformedfel.GetGramShape();
    const double absdet = transformedfel.GetAbsDeterminant();
    const double d = local_coefs.GetDiffusionCoef(pPart)/local_coefs.GetHenryWeighting(pPart);
    const double hw = 1./local_coefs.GetHenryWeighting(pPart);
	
    const Quad3CL<Point3DCL> & q3_u = local_coefs.GetVelocityAsQuad3();

    if (transformedfel.stabilized()){
      StabilizedTransformedP1FiniteElement& stabfe = static_cast<StabilizedTransformedP1FiniteElement&>(transformedfel);      
      for(int i= 0; i < 4; ++i) {
          Quad3CL<> & phi_i = stabfe.GetBaseShapeAsQuad3CL(i);
          Quad3CL<> & beta_i = stabfe.GetStabTestShapeAsQuad3CL(i);
          Quad3CL<> & phi_plus_beta_i = stabfe.GetTestShapeAsQuad3CL(i);
          for(int j= 0; j < 4; ++j) {
              Quad3CL<> & phi_j = stabfe.GetBaseShapeAsQuad3CL(j);
              Quad3CL<> & beta_j = stabfe.GetStabTestShapeAsQuad3CL(j);
              elmats.M[i][j] = hw*Quad3CL<>( phi_plus_beta_i*phi_j).quad( absdet);
              //elmats.Mr[i][j]= hw*Quad3CL<>( phi_plus_beta_i*phi_j*local_coefs.GetReactionAsQuad3()).quad( absdet);
              elmats.A[i][j] = d*GTG( i, j)/6.0*absdet;
              elmats.A[i][j]+= hw*Quad3CL<>( beta_i*beta_j).quad( absdet)*stabfe.GetStabilizationParameter();
              elmats.C[i][j] = hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( j)))*phi_i).quad( absdet);
          }
      }
    }
    else
    {  
      for(int i= 0; i < 4; ++i) {
          Quad3CL<> & phi_i = transformedfel.GetGridfunctions().GetShapeAsQuad3(i);
          for(int j= 0; j < i; ++j) {
              Quad3CL<> & phi_j = transformedfel.GetGridfunctions().GetShapeAsQuad3(j);
              elmats.M[j][i]= hw*P1DiscCL::GetMass( i, j)*absdet;
              elmats.M[i][j]= elmats.M[j][i];
              elmats.A[j][i]= d*GTG( i, j)/6.0*absdet;
              elmats.A[i][j]= elmats.A[j][i];
              elmats.C[i][j]= hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( j)))*phi_i).quad( absdet);
              elmats.C[j][i]= hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( i)))*phi_j).quad( absdet);
          }
          elmats.M[i][i]= hw*P1DiscCL::GetMass( i, i)*absdet;
          elmats.A[i][i]= d*GTG( i, i)/6.0*absdet;
          elmats.C[i][i]= hw*Quad3CL<>( dot( q3_u, Quad3CL<Point3DCL>( G.col( i)))*phi_i).quad( absdet);
      }
    }

}

// Couplings between the XFEM basis functions (at the same time step or different time steps) when the tetra is cut by only one interface.
void SetupLocalOneInterfaceSystem( TransformedP1FiniteElement& transformedfel, InterfaceTetraCL& cut, 
      ConvDiffElementMatrices & elmats, LocalConvDiffReacCoefficients& local_coefs)
{
    const SMatrixCL<3,4> G = transformedfel.GetDShape();
    const SMatrixCL<4,4> GTG = transformedfel.GetGramShape();
    const double absdet = transformedfel.GetAbsDeterminant();

    StabilizedTransformedP1FiniteElement* stabfe = NULL;
    if (transformedfel.stabilized())
      stabfe = static_cast<StabilizedTransformedP1FiniteElement*>(&transformedfel);
      
    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumTetra(); //# of subtetras

    elmats.ResetSigned(); //CL: <--- kann weg
    
    for (Uint k=0; k< NumTets; ++k){
        bool IAmInPosPart = k>=cut.GetNumNegTetra();   
		
        double d = local_coefs.GetDiffusionCoef(IAmInPosPart)/local_coefs.GetHenryWeighting(IAmInPosPart);
        double hw = 1./local_coefs.GetHenryWeighting(IAmInPosPart);
           
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        
        if (!IsRegBaryCoord(T))  continue;
        
        transformedfel.SetSubTetra(T);
        if (stabfe) stabfe->CalcStabilization(IAmInPosPart);
        
        double Vol = absdet*VolFrac(T);
        if (isnan(Vol)|| isinf(Vol)){
            std::cout<<"Vol " <<VolFrac(T)<<"\n";
            std::cout<<"SetupLocalOneInterfaceSystem: Support of XFEM is too small.\t";
            continue;
        }
        Quad3CL<Point3DCL> q3_u(local_coefs.GetVelocityAsLocalP2(), transformedfel.GetNodes());
        bool irreg = false; 

        for(int i= 0; i < 4; ++i) {
            Quad3CL<> qp1(transformedfel.GetGridfunctions().GetShapeAsLocalP1(i), transformedfel.GetNodes());
            for(int j= 0; j < 4; ++j) {
              
                double iM = 0;
                
                Quad3CL<> qC(dot( q3_u, Quad3CL<Point3DCL>(G.col( j)))*qp1); // G is constant in t
                double iC = qC.quad(Vol) * hw; 
                double iA = Vol* GTG( j, i)/6. * d; 
                
                if (stabfe){
                  Quad3CL<> qM(stabfe->GetBaseShapeAsQuad3CL(j)*stabfe->GetTestShapeAsQuad3CL(i));
                  iM = qM.quad(Vol) * hw;
                  Quad3CL<> qSC(stabfe->GetStabTestShapeAsQuad3CL(j)*stabfe->GetStabTestShapeAsQuad3CL(i));
                  iA += qSC.quad(Vol) * hw * stabfe->GetStabilizationParameter();
                }else
                {
                  Quad3CL<> qM(transformedfel.GetGridfunctions().GetProductShapeAsLocalP2(i,j),transformedfel.GetNodes());
                  iM = qM.quad(Vol) * hw;
                }              
              
                if (isnan(iM)|| isinf(iM)||isnan(iA)|| isinf(iA)||isnan(iC)|| isinf(iC)) {
                    elmats.ResetSigned();
                    irreg = true;
                    break;
                }
       
              
       
                // D and H are piecewise constants in T. When
                // - Compute the integrals between XFEM basis functions at the same time step
                // - Compute the integrals between new FEM basis functions with old XFEM
                // basis functions (wrt old interface) -> (FEM function, old XFEM).
       
                if (!IAmInPosPart){
                    elmats.M_n[i][j]+= iM;
                    elmats.A_n[i][j]+= iA;
                    elmats.C_n[i][j]+= iC;
                }
                else{
                    elmats.M_p[i][j]+= iM;
                    elmats.A_p[i][j]+= iA;
                    elmats.C_p[i][j]+= iC;
                }
            }
            if (irreg) break;
        }
    }
}

void SetupLocalOnePhaseMassMatrix( double locM[4][4], TransformedP1FiniteElement & transfp1fel, const double H, bool pPart)
{
    const double h = pPart ? 1. : 1./H;
    
    if (transfp1fel.stabilized()){
      StabilizedTransformedP1FiniteElement* stabfe = static_cast<StabilizedTransformedP1FiniteElement*>(&transfp1fel);
      for(int i= 0; i < 4; ++i) {
        for(int j= 0; j < 4; ++j) {
            locM[i][j] = h * Quad3CL<>(stabfe->GetBaseShapeAsQuad3CL(j)*stabfe->GetTestShapeAsQuad3CL(i)).quad(transfp1fel.GetAbsDeterminant());
        }
      }
    }
      else
    {
      for(int i= 0; i < 4; ++i) {
          for(int j= 0; j < i; ++j) {
              locM[j][i]= h*P1DiscCL::GetMass( i, j)*transfp1fel.GetAbsDeterminant();
              locM[i][j]= locM[j][i];
          }
          locM[i][i]= h*P1DiscCL::GetMass( i, i)*transfp1fel.GetAbsDeterminant();
      }
    }
}

/// compute the mass matrix for the case that the tetrahedron is cut ONLY for the old OR new time. \n
/// computes the off-diagonal block of \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ M_{1,2}[i,j] = \int h \phi_i^{new} \phi_j^{old} \f$ \n
/// or
/// \f$ M_{2,1}[i,j] = \int h \phi_i^{old} \phi_j^{new} \f$ \n
void SetupLocalOneInterfaceMassMatrix( InterfaceTetraCL& cut, double M_n[4][4], double M_p[4][4], 
    TransformedP1FiniteElement & transfp1fel, const double H, bool sign[4], bool jumps, bool pPart)
{
    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumTetra(); /// # of subtetras
    // D, H are constant if T \cap Gamma_old is empty
    double h = pPart ? 1. : 1./H;
    std::memset( M_n,0, 4*4*sizeof(double));
    std::memset( M_p,0, 4*4*sizeof(double));
    
    StabilizedTransformedP1FiniteElement* stabfe = NULL;
    if (transfp1fel.stabilized())
      stabfe = static_cast<StabilizedTransformedP1FiniteElement*>(&transfp1fel);
    
    for(int i= 0; i < 4; ++i) {
        sign[i]= (cut.GetSign(i) == 1);  
    }        

    for (Uint k=0; k< NumTets; ++k){
        const SArrayCL<BaryCoordCL,4>& T =cut.GetTetra(k);
        transfp1fel.SetSubTetra(T);
        if (stabfe) stabfe->CalcStabilization(k>=cut.GetNumNegTetra());
        if (!IsRegBaryCoord(T)) continue;
        
        double Vol = transfp1fel.GetAbsDeterminant()*VolFrac(T);
        if (isnan(Vol)|| isinf(Vol)){
            std::cout<<"Vol " <<VolFrac(T)<<"\n";
            std::cout<<" Support of XFEM is too small.\t";
            continue;
        }
        bool irreg=false;
        for(int i= 0; i < 4; ++i) {
            for(int j= 0; j < 4; ++j) {
                double iM = 0;
                if (stabfe){
                  Quad3CL<> qM(stabfe->GetBaseShapeAsQuad3CL(j)*stabfe->GetTestShapeAsQuad3CL(i));
                  iM = qM.quad(Vol);
                }else
                {
                  Quad3CL<> qM(transfp1fel.GetBaseShapeAsQuad3CL(j)*transfp1fel.GetBaseShapeAsQuad3CL(i));
                  //Quad3CL<> qM(transfp1fel.GetGridfunctions().GetProductShapeAsLocalP2(i,j),transfp1fel.GetNodes());
                  iM = qM.quad(Vol);
                }
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
                        M_n[i][j]+= iM/H;
                    else
                        M_p[i][j]+= iM ;
                }
                // D and H are constant.
                // Compute (new XFEM test function, FEM) at old time step. Tetra is cut by only new interface.
                else{  
                    if (k<cut.GetNumNegTetra())
                        M_n[i][j]+= iM*h;
                    else
                        M_p[i][j]+= iM*h;
                }
            }
            if (irreg) break;  //leave the i-loop
        }
    }
}

/// compute the mass matrix for the case that the tetrahedron is cut for the old and new time. \n
/// computes \f$ M(u,v) = \int h u v \f$ with h inverse of the henry coefficient of the domain \n
/// \f$ M_{2,1}[i,j] = \int h \phi_i^{old} \phi_j^{new} \f$ \n
/// \f$ M_{2,2}[i,j] = \int h \phi_i^{old} \phi_j^{old} \f$ \n
void SetupLocalTwoInterfacesMassMatrix( InterfaceTetraCL& cut, InterfaceTetraCL& oldcut, 
    double M22[4][4], double M21[4][4], TransformedP1FiniteElement & transfp1fel, const double H, LocalP2CL<>& lp2_oldlset)
{
    bool sign[4], oldsign[4];
    std::memset( M22, 0, 4*4*sizeof(double));
    std::memset( M21, 0, 4*4*sizeof(double));
    
    for(int i= 0; i < 4; ++i) {
        sign[i]= (cut.GetSign(i) == 1);
        oldsign[i]= (oldcut.GetSign(i) == 1);
    }        

    StabilizedTransformedP1FiniteElement* stabfe = NULL;
    if (transfp1fel.stabilized())
      stabfe = static_cast<StabilizedTransformedP1FiniteElement*>(&transfp1fel);

    cut.ComputeSubTets();
    Uint NumTets=cut.GetNumTetra(); /// # of subtetras

    for (Uint k=0; k< NumTets; ++k){
        bool Tk= (k>=cut.GetNumNegTetra());    // Tk in Omega_new_+?
        const SArrayCL<BaryCoordCL,4>& T =  cut.GetTetra(k);
        if (!IsRegBaryCoord(T)) continue;
        
        transfp1fel.SetSubTetra(T);
        if (stabfe) stabfe->CalcStabilization(Tk);
        
        InterfaceTetraCL suboldcut;
        suboldcut.Init(T, lp2_oldlset,0.);
        double VolT = transfp1fel.GetAbsDeterminant()*VolFrac(T);
        bool nocut= !suboldcut.Intersects();
        suboldcut.ComputeSubTets();
        Uint NumOldTets= suboldcut.GetNumTetra();
        if (nocut){  
            bool pPart = (suboldcut.GetSign( 0) == 1);
            double h = pPart ? 1. : 1./H;

            bool irreg=false;
            for(int i= 0; i < 4; ++i) {
                if (Tk == sign[i]) continue; // supp(Phi_i^{Gamma_new}) \cap Tk is empty
                for(int j= 0; j < 4; ++j) {
                    double iM = 0;
                    if (stabfe){
                      Quad3CL<> qM(stabfe->GetBaseShapeAsQuad3CL(j)*stabfe->GetTestShapeAsQuad3CL(i));
                      iM = qM.quad(VolT)*h;
                    }
                    else
                    {
                      Quad3CL<> qM(transfp1fel.GetBaseShapeAsQuad3CL(i)*transfp1fel.GetBaseShapeAsQuad3CL(j));
                    //Quad3CL<> qM(transfp1fel.GetGridfunctions().GetProductShapeAsLocalP2(i,j),transfp1fel.GetNodes());
                      iM = qM.quad(VolT)*h;
                    }                  
                  
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
                    M21[i][j]+= sign[i]? -iM: iM;

                    if (pPart== oldsign[j]) continue; // supp(Phi_j^{Gamma_old}) \cap Tk is empty

                    M22[i][j]+= (sign[i]==oldsign[j]) ? iM : -iM ;
                }
                if (irreg) break;
            }
            continue;
        }
           
        for (Uint m=0; m< NumOldTets; ++m){
            bool Tkm= (m>=suboldcut.GetNumNegTetra());  // Tkm in Omega_old_+?
            double h = Tkm ? 1. : 1./H;
            const SArrayCL<BaryCoordCL,4>& Tc =  suboldcut.GetTetra(m);
            if (!IsRegBaryCoord(Tc)) continue;
            
            transfp1fel.SetSubTetra(Tc);
            if (stabfe) stabfe->CalcStabilization(Tkm);            
            
            double Vol = transfp1fel.GetAbsDeterminant()*VolFrac(Tc);
            bool irreg=false;
            for(int i= 0; i < 4; ++i) {
                if (Tk == sign[i]) continue; // supp(Phi_i^{Gamma_new}) \cap Tk is empty
                for(int j= 0; j < 4; ++j) {
                  
                    double iM = 0;
                    if (stabfe){
                      Quad3CL<> qM(stabfe->GetBaseShapeAsQuad3CL(j)*stabfe->GetTestShapeAsQuad3CL(i));
                      iM = qM.quad(Vol)*h;
                    }else
                    {
                      Quad3CL<> qM(transfp1fel.GetBaseShapeAsQuad3CL(i)*transfp1fel.GetBaseShapeAsQuad3CL(j));
//                    Quad3CL<> qM(transfp1fel.GetGridfunctions().GetProductShapeAsLocalP2(i,j),transfp1fel.GetNodes());
                      iM = qM.quad(Vol)*h;
                    }                          
                  
                    if (isnan(iM)|| isinf(iM)) {
                        ///> TODO: If a local value in a child tetrahedron is irregular, ignore this tetra
                        ///> The contribution of other tetra must be preserved
                        ///> For example tmpM21[4][4]
//                         std::memset( M21,0, 4*4*sizeof(double));
//                         std::memset( M22,0, 4*4*sizeof(double));
                        irreg=true;
                        break; //leave the j-loop 
                    }
                    M21[i][j]+= sign[i]? -iM: iM;
                    if (Tkm == oldsign[j]) continue;// supp(Phi_j^{Gamma_old}) \cap Tkm is empty

                    M22[i][j]+= (sign[i]==oldsign[j]) ? iM : -iM ;
                }
                if (irreg) break;
            }
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

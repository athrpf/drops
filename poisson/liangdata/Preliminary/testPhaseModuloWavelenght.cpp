#include <iostream>
#include "boxinterp.hpp"

int main(){

  int NX = 366;
  double deltaX = 0.00002;
  double L = NX*deltaX;
  double phi0=-1;
  double phi;
  double phiModuloWavelenght;
  double errormax=1.e-14;
  
  while(phi0 <= 0 || phi0 >= L){
  
  std::cout << "Enter a phase-value phi0, that is >" << 0 << " and <" <<  L << std::endl;
  
  std::cin >> phi0;
  
  }
  
  std::cout << " phiModuloWavelenght - phi0 should be zero. An error report will occur, if abs(phiModuloWavelenght - phi0) is larger than " << errormax  << std::endl;
  
  for(int k=0; k<10000; k++){
  
               phi = phi0 + k*L;
  
               phiModuloWavelenght = phaseModuloWavelenght(phi, L);
          
       if(phiModuloWavelenght - phi0 < -errormax || phiModuloWavelenght - phi0 > errormax){
          std::cout << "Error for k=" << k << " :" << "phiModuloWavelenght - phi0 =" << phiModuloWavelenght - phi0 << ", not equal to zero."     
          << std::endl;
       }
  }
  
  return 0;
}

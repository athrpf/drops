/// \file interfacepatch_doublecut.cpp
/// \brief tests implementation of the numbering methods of interface finite elements
/// \author LNM RWTH Aachen: Christoph Lehrenfeld

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

#include "geom/simplex.h"
#include "levelset/levelset.h"
#include "num/interfacePatch.h"
#include <fstream>

using namespace DROPS;

double plane(const Point3DCL& p)
{
    return p[0] - 0.4;
}

double plane2(const Point3DCL& p)
{
    return p[0] - 0.7;
}

double plane3(const Point3DCL& p)
{
    return p[0]*p[1] - 0.05;
}

double sphere(const Point3DCL& p)
{
    return p.norm_sq()-1.0/16.0;
}


void TestVolumeOnRefTetra()
{
    //make the reference tetraeder
    Point3DCL p0(0.), p1(0.), p2(0.), p3(0.);
    p1[0]= p2[1]= p3[2]= 1;
    TetraBuilderCL tet (0,p0, p1, p2, p3);
    MultiGridCL MG( tet);
    TetraCL & reftet = *MG.GetTetrasBegin();

    LocalP2CL<double> ls1(0.);
    LocalP2CL<double> ls2(0.);
    for (int i=0; i<10; i++){
        ls1[i]= sphere(i<4 ? reftet.GetVertex(i)->GetCoord() : GetBaryCenter( *reftet.GetEdge(i-4)));
        ls2[i]= plane2(i<4 ? reftet.GetVertex(i)->GetCoord() : GetBaryCenter( *reftet.GetEdge(i-4)));
    }

    InterfaceTetraCL firstcut;
    firstcut.Init( reftet, ls1, 0.0);

    firstcut.ComputeSubTets();

    Uint firstNumTets=firstcut.GetNumTetra(); /// # of subtetras

    double controlvol1 = 0.0;
    double controlvol2 = 0.0;
    double controlvol1_p = 0.0;
    double controlvol2_p = 0.0;
    double controlvol1_m = 0.0;
    double controlvol2_m = 0.0;
    
    for (Uint k=0; k< firstNumTets; ++k){
        bool pPart_first = (k>=firstcut.GetNumNegTetra());

        const SArrayCL<BaryCoordCL,4>& subtet =  firstcut.GetTetra(k);
        controlvol1 += VolFrac(subtet);
        if (pPart_first)
            controlvol1_p += VolFrac(subtet);
        else
            controlvol1_m += VolFrac(subtet);

        InterfaceTetraCL secondcut;

        Uint child = firstcut.GetChildIdx(k);
        secondcut.Init( subtet, ProjectIsoP2ChildToParentP1(ls2,child), 0.0);
        secondcut.ComputeSubTets(/* subdivide_first = */false);
        double subvol = 0.0;
        double subtetvol = VolFrac(subtet);
        if (secondcut.Intersects()){
            Uint secondNumTets=secondcut.GetNumTetra(); /// # of subtetras
            for (Uint j=0; j< secondNumTets; ++j){
                bool pPart_second = (j>=secondcut.GetNumNegTetra());
                const SArrayCL<BaryCoordCL,4>& subsubtet = secondcut.GetTetra(j);
                subvol += VolFrac(subsubtet);
                if (pPart_second)
                    controlvol2_p += VolFrac(subsubtet);
                else
                    controlvol2_m += VolFrac(subsubtet);

            }
            if ((std::fabs(subvol-subtetvol)/subtetvol)>1e-8 ){
                std::cout << " PROBLEMS on subtet " << k << std::endl;
                std::cout << " ACCUMULATED subsubtet volume is " << subvol/6.0 << std::endl;
                std::cout << " while subtet volume is " << subtetvol/6.0 << std::endl;
            }
        }
        else{
            subvol = subtetvol;
            bool pPart_second  = (secondcut.GetSign( 0) == 1);
            
            if (pPart_second)
                controlvol2_p += VolFrac(subtet);
            else
                controlvol2_m += VolFrac(subtet);
        }
        controlvol2 += subvol;
    }  

    std::cout << "ACCUMULATED    subtet volume = " << controlvol1 << " (should be 1.0) " << std::endl;
    std::cout << "ACCUMULATED subsubtet volume = " << controlvol2 << " (should be 1.0) " << std::endl;

    std::cout << "ACCUMULATED    subtet volume (plus ) = " << controlvol1_p << " (should be ...) " << std::endl;
    std::cout << "ACCUMULATED subsubtet volume (plus ) = " << controlvol2_p << " (should be ...) " << std::endl;
    std::cout << "ACCUMULATED    subtet volume (minus) = " << controlvol1_m << " (should be ...) " << std::endl;
    std::cout << "ACCUMULATED subsubtet volume (minus) = " << controlvol2_m << " (should be ...) " << std::endl;

}


int main ()
{
  try 
    {
      TestVolumeOnRefTetra();
      std::cout << "Tested Volume on RefTetra" << std::endl;
      return 0;
    }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

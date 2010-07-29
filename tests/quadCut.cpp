/// \file quadCut.cpp
/// \brief tests quadrature on cuts
/// \author LNM RWTH Aachen: Sven Gross, Volker Reichelt; SC RWTH Aachen: Oliver Fortmeier

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

#include "geom/builder.h"
#include "levelset/levelset.h"

using namespace DROPS;

void ComputeIntOnCuts( const TetraCL& t, const VecDescCL& ls, const LocalP2CL<>& f, double absdet, VectorCL& quadPos, VectorCL& quadNeg)
{
    InterfaceTetraCL cut;

    BndCondT bc[6]= { NoBC, NoBC, NoBC, NoBC, NoBC, NoBC };
    BndDataCL<>::bnd_val_fun bfun[6]= {0,0,0,0,0,0};
    BndDataCL<> bnd(6, bc, bfun);

    cut.Init( t, ls, bnd);
    if (!cut.Intersects()) // kein Beitrag
    { std::cout << "No cut on tetra...\n"; }

    for (int ch=0; ch<8; ch++)
    {
        cut.ComputeCutForChild(ch);
        if (ch==0) cut.DebugInfo( std::cout, true);
        quadPos[ch]= cut.quad( f, absdet, true);  // integrate on positive part
        quadNeg[ch]= cut.quad( f, absdet, false); // integrate on negative part
    }
}

bool CheckSum( const VectorCL& v1, const VectorCL& v2, const VectorCL& v)
{
    // check v == v1 + v2
    const VectorCL diff(v1+v2-v);
    if (norm(diff)>1e-18)
    {
        std::cout << ">>> inconsistency!\nquadPos = \t" << v1 << "\nquadNeg = \t" << v2 << "\nsum = \t\t" << VectorCL(v1+v2)
                  << "should be\t" << v << "\ndiff = \t\t" << diff << std::endl;
        return false;
    }
    return true;
}

int main ()
{
  try {
    std::cout << "Testing quadrature on cuts of tetra." << std::endl;
    TetraBuilderCL tetra( 0); // unrefined reference tetra
    MultiGridCL mg( tetra);
    TetraCL& t= *mg.GetAllTetraBegin();

    IdxDescCL idx( P2_FE);
    idx.CreateNumbering( 0, mg);
    VecDescCL ones( &idx), ls( &idx);
    for (int i=0; i<10; ++i)
        ones.Data[i]= 1;

    const double absdet= 1;
    bool ok= true;
    VectorCL quadSum( 8), quadPos( 8), quadNeg( 8); // integrals on 8 children

    for (int i=0; i<10; ++i)
    {
        LocalP2CL<> f;    f[i]= 1.;
        std::cout << "======== testing P2 basis function " << i << "===========\n";

        ComputeIntOnCuts( t, ones, f, absdet, quadSum, quadNeg);
        if (norm(quadNeg)>0) { std::cout << ">>> quadNeg should be zero!\n"; return 1; }

        Quad2CL<> q( f);
        const double integral= q.quad(absdet),
            diff= std::abs(integral-quadSum.sum());
        if (diff>1e-17)
        {
            std::cout << ">>> sum of quadSum should be equal to integral over tetra!\n";
            std::cout << quadSum.sum() << " != " << integral << ",\tdiff = " << diff << std::endl;
            return 1;
        }

        ls.Data= ones.Data; ls.Data[0]= -3;
        ComputeIntOnCuts( t, ls, f, absdet, quadPos, quadNeg);
        ok= ok && CheckSum( quadPos, quadNeg, quadSum);

        ls.Data[4]= -7;
        ComputeIntOnCuts( t, ls, f, absdet, quadPos, quadNeg);
        ok= ok && CheckSum( quadPos, quadNeg, quadSum);

        ls.Data[5]= 0;
        ComputeIntOnCuts( t, ls, f, absdet, quadPos, quadNeg);
        ok= ok && CheckSum( quadPos, quadNeg, quadSum);

        ls.Data[0]= 0;
        ComputeIntOnCuts( t, ls, f, absdet, quadPos, quadNeg);
        ok= ok && CheckSum( quadPos, quadNeg, quadSum);

        ls.Data[7]= 0;
        ComputeIntOnCuts( t, ls, f, absdet, quadPos, quadNeg);
        ok= ok && CheckSum( quadPos, quadNeg, quadSum);
    }
    if (ok) std::cout << "\n\nCONGRATULATIONS! All tests successfully!\n";
    else    std::cout << "\n\n>>> HELP! Some errors occured!\n";
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

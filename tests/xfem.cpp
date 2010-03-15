/// \file xfem.cpp
/// \brief tests XFEM-functions with a planar interface
/// \author LNM RWTH Aachen: Patrick Esser, Joerg Grande, Sven Gross, Trung Hieu Nguyen; SC RWTH Aachen: Oliver Fortmeier

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

#include "misc/utils.h"
#include "num/spmat.h"
#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "num/discretize.h"
#include "num/fe.h"
#include "misc/problem.h"
#include "levelset/levelset.h"
#include "levelset/surfacetension.h"
#include "stokes/instatstokes2phase.h"
#include "out/output.h"
#include "out/ensightOut.h"
#include <fstream>
#include <iomanip>

using namespace DROPS;

inline int iSign( double x)
{
    return x<0 ? -1 : x>0 ? 1 : 0;
}

// \Omega_1 is the domain with phasebnd < 0.
double phasebnd (const Point3DCL& p)
{
    return p[1] + p[2] - 0.05;
//    return p[1] + p[2] - 0.55;
//    return p[0] + p[1] - 0.5;
}

double f1 (const Point3DCL& p, double)
{
    return p.norm_sq();
}

double f2 (const Point3DCL& p, double)
{
    return 3.0*p[0]*p[0] + p[1]*p[1] + 2.0*p[2]*p[2] + 2.0;
}

// double f1 (const Point3DCL& p, double)
// {
//     return p[0] + 0.5*p[1]  - 1.0;
// }
//
// double f2 (const Point3DCL& p, double)
// {
//     return 0.7*p[0] + p[1] + 0.3*p[2];
// }

double f1sq (const Point3DCL& p, double)
{
    return std::pow( f1( p, 0.), 2);
}

double f2sq (const Point3DCL& p, double)
{
    return std::pow( f2( p, 0.), 2);
}

void InitPiecewiseP2 (instat_scalar_fun_ptr fneg, VecDescCL& un,
    instat_scalar_fun_ptr fpos, VecDescCL& up, MultiGridCL& MG, VecDescCL& Phi)
{
    const Uint idx= Phi.RowIdx->GetIdx(),
               unidx= un.RowIdx->GetIdx(),
               upidx= up.RowIdx->GetIdx();


    DROPS_FOR_TRIANG_VERTEX( MG, MG.GetLastLevel(), it) {
        switch (iSign( Phi.Data[it->Unknowns(idx)])) {
          case 1:
            un.Data[it->Unknowns( unidx)]= 0.;
            up.Data[it->Unknowns( upidx)]= fpos( it->GetCoord(), 0.);
            break;
          case -1:
            un.Data[it->Unknowns( unidx)]= fneg( it->GetCoord(), 0.);
            up.Data[it->Unknowns( upidx)]= 0.;
            break;
          default:
            un.Data[it->Unknowns( unidx)]= fneg( it->GetCoord(), 0.);
            up.Data[it->Unknowns( upidx)]= fpos( it->GetCoord(), 0.);
        }
    }
    DROPS_FOR_TRIANG_EDGE( MG, MG.GetLastLevel(), it) {
        switch (iSign( Phi.Data[it->Unknowns(idx)])) {
          case 1:
            un.Data[it->Unknowns( unidx)]= 0.;
            up.Data[it->Unknowns( upidx)]= fpos( GetBaryCenter( *it), 0.);
            break;
          case -1:
            un.Data[it->Unknowns( unidx)]= fneg( GetBaryCenter( *it), 0.);
            up.Data[it->Unknowns( upidx)]= 0.;
            break;
          default:
            un.Data[it->Unknowns( unidx)]= fneg( GetBaryCenter( *it), 0.);
            up.Data[it->Unknowns( upidx)]= fpos( GetBaryCenter( *it), 0.);
        }
    }
    InterfacePatchCL cut;
    DROPS_FOR_TRIANG_TETRA( MG, MG.GetLastLevel(), it) {
        cut.Init( *it, Phi);
        if (cut.Intersects()) {
            for (Uint i= 0; i < 4; ++i) {
                const VertexCL& v= *it->GetVertex( i);
                un.Data[v.Unknowns( unidx)]= fneg( v.GetCoord(), 0.);
                up.Data[v.Unknowns( upidx)]= fpos( v.GetCoord(), 0.);
            }
            for (Uint i= 0; i < 6; ++i) {
                const EdgeCL& e= *it->GetEdge( i);
                un.Data[e.Unknowns( unidx)]= fneg( GetBaryCenter( e), 0.);
                up.Data[e.Unknowns( upidx)]= fpos( GetBaryCenter( e), 0.);
            }
        }
    }
}

void P1XOnPart (const VecDescCL& p1x, const ExtIdxDescCL& Xidx, VecDescCL& p_part,
    const LevelsetP2CL& lset, bool posPart)
{
    const Uint lvl= p1x.RowIdx->TriangLevel(),
        idxnum= p1x.RowIdx->GetIdx();
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    const MultiGridCL& mg= lset.GetMG();

    p_part.SetIdx( p1x.RowIdx);
    VectorCL& pp= p_part.Data;
    pp= p1x.Data; // Our Ensight-writer ignores the extra p1x values at the end.

    // add extended pressure
    DROPS_FOR_TRIANG_CONST_VERTEX( mg, lvl, it) {
        const IdxT nr= it->Unknowns( idxnum);
        if (Xidx[nr] == NoIdx) continue;

        const bool is_pos= InterfacePatchCL::Sign( ls.val( *it)) == 1;
        if (posPart == is_pos) continue; // extended hat function == 0 on this part
        if (posPart)
            pp[nr]+= p1x.Data[Xidx[nr]];
        else
            pp[nr]-= p1x.Data[Xidx[nr]];
    }
}

class CoeffCL
{
  public:
    double mu (double) const { return 1.0; }
};

int main (int argc, char** argv)
{
  try {
    int numref=10;
    double xfemstab=0.;
    if (argc==3) {
        numref=atoi(argv[1]);
        xfemstab=atof(argv[2]);
    }

    std::cout << "numref: " << numref << "\txfemstab: " << xfemstab <<'\n';

    DROPS::BrickBuilderCL brick( Point3DCL( -1.0),
                                  2.0*std_basis<3>(1),
                                  2.0*std_basis<3>(2),
                                  2.0*std_basis<3>(3),
                                  numref, numref, numref);
    MultiGridCL mg( brick);

//    TetraBuilderCL builder( 0);//,   std_basis<3>( 0), 2*std_basis<3>( 1),
                               //2*std_basis<3>( 2), 2*std_basis<3>( 3));
 //   MultiGridCL mg( builder);

    instat_scalar_fun_ptr sigma (0);
    SurfaceTensionCL sf( sigma, 0);
    LevelsetP2CL lset( mg, sf);
    lset.idx.CreateNumbering( mg.GetLastLevel(), mg);
    lset.Phi.SetIdx( &lset.idx);
    lset.Init( &phasebnd);

    IdxDescCL idx( P1X_FE, BndCondCL(0), 0, /*omit_bound=*/ xfemstab);
    idx.CreateNumbering( mg.GetLastLevel(), mg, &lset.Phi);
    ExtIdxDescCL& extidx= idx.GetXidx();
    std::cout << "P1-Unbekannte: " << extidx.GetNumUnknownsStdFE()
              << " P1X-Unbekannte: " << idx.NumUnknowns() << '\n';

    VecDescCL beta( &idx), b( &idx);

    // For the ensight visualisation of the piecewise quadratic function f1-f2
    IdxDescCL p2idx( P2_FE);
    p2idx.CreateNumbering( mg.GetLastLevel(), mg);
    VecDescCL uneg( &p2idx), upos( &p2idx);
    InitPiecewiseP2( f1, uneg, f2, upos, mg, lset.Phi);

    // Setup the mass matrix
    MatDescCL M( &idx, &idx);
    SetupPrMass_P1X( mg, CoeffCL(), M.Data, idx, lset);

    // Setup the right hand side
    IdxT Numb[4];
    double absdet;
    InterfaceTetraCL cut;
    Quad5CL<> qf;

    DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), sit) {
        GetLocalNumbP1NoBnd( Numb, *sit, idx);
        absdet= sit->GetVolume()*6.0;
        cut.Init( *sit, lset.Phi);

        cut.ComputeSubTets();
        BaryCoordCL* nodes;
        for (Uint i= 0; i < 4; ++i) {
            const IdxT xidx= extidx[Numb[i]];
            const bool have_xidx( xidx != NoIdx);
            double intp1= 0.0, intp1x= 0.0, intpos= 0.0, intneg= 0.0;
            LocalP1CL<> phip1;
            phip1[i]= 1.0;
            for (Uint k=0; k< cut.GetNumTetra(); ++k)
            {
                nodes = qf.TransformNodes(cut.GetTetra(k));
                qf.assign(*sit, cut.GetNumNegTetra()>k ? &f1 : &f2, 0.0, nodes);
                Quad5CL<> qphi(phip1, nodes);
                if (cut.GetNumNegTetra()>k)
                    intneg += Quad5CL<>(qphi*qf).quad(absdet*VolFrac(cut.GetTetra(k)));
                else
                    intpos += Quad5CL<>(qphi*qf).quad(absdet*VolFrac(cut.GetTetra(k)));
                delete[] nodes;
            }
            intp1+= intpos + intneg;
            if (have_xidx)
                intp1x+= cut.GetSign( i) == 1 ? -intneg : intpos;
            b.Data[Numb[i]]+= intp1;
            if (have_xidx) b.Data[xidx]+= intp1x;
        }
    }

    double intval=0.;
    DROPS_FOR_TRIANG_TETRA( mg, mg.GetLastLevel(), sit) {
        absdet= sit->GetVolume()*6.0;
        cut.Init( *sit, lset.Phi);
        cut.ComputeSubTets();
        BaryCoordCL* nodes;
        for (Uint i=0; i<cut.GetNumTetra(); i++)
        {
            nodes = qf.TransformNodes(cut.GetTetra(i));
            qf.assign( *sit, cut.GetNumNegTetra()>i ? &f1sq : &f2sq, 0.0, nodes);
            intval += qf.quad(absdet*std::fabs( VolFrac(cut.GetTetra(i))));
            delete[] nodes;
        }
    }

    // Solve the linear system...
    int max_iter= 200;
    double tol= 1e-16;
    PCG( M.Data, beta.Data, b.Data, JACPcCL( 1.0), max_iter, tol, /*measure_relative_tol*/ true);
    std::cout <<  "iter: " << max_iter << "\ttol: " << tol << '\n';

    //Ensight output
    NoBndDataCL<> ubnd;
    Ensight6OutCL ensight( "xfem.case", 0);
    const std::string filename= "ensight/xfem";
    ensight.Register( make_Ensight6Geom     ( mg, mg.GetLastLevel(),           "Cube",     filename + ".geo"));
    ensight.Register( make_Ensight6Scalar   ( lset.GetSolution(),              "Levelset", filename + ".scl"));
    ensight.Register( make_Ensight6Scalar   ( make_P2Eval( mg, ubnd, upos),    "up",       filename + ".up"));
    ensight.Register( make_Ensight6Scalar   ( make_P2Eval( mg, ubnd, uneg),    "un",       filename + ".un"));
    // Output the L_2-projection
    ensight.Register( make_Ensight6P1XScalar( mg, lset.Phi, beta,              "ul",       filename + ".ul"));

    ensight.Write();

    //1D-Plots
    VecDescCL ulneg, ulpos;
    P1XOnPart( beta, extidx, ulpos, lset, true);
    P1XOnPart( beta, extidx, ulneg, lset, false);
    std::ofstream out ("u.txt");
    std::ofstream outpr ("up.txt");
    Point3DCL p;
    LevelsetP2CL::const_DiscSolCL ls= lset.GetSolution();
    DROPS_FOR_TRIANG_VERTEX( mg, mg.GetLastLevel(), it) {
        p= it->GetCoord();
        if (p[0]==0 && p[1]==0) {
            const double u=  lset.Phi.Data[it->Unknowns( lset.idx.GetIdx())] > 0. ? upos.Data[it->Unknowns( p2idx.GetIdx())] : uneg.Data[it->Unknowns( p2idx.GetIdx())];
            const IdxT nr= it->Unknowns( idx.GetIdx());
            const double up= lset.Phi.Data[it->Unknowns( lset.idx.GetIdx())] > 0. ? ulpos.Data[nr] : ulneg.Data[nr];
            out << p[2] << " " << u << '\n';
            outpr << p[2] << " " << up << '\n';
        }
    }

    std::cout << std::setprecision(20);
    std::cout << "||u_l||_0^2 = " << dot (M.Data*beta.Data, beta.Data) <<'\n';
    std::cout << "|| u ||_0^2 = " << intval << '\n';
    const double err = intval - dot (M.Data*beta.Data, beta.Data);
    std::cout << "|| u - u_l ||_0 = " << std::scientific << std::sqrt(std::abs(err)) << '\n';
    if (err<0) std::cout << "Fehler: Norm der Projektion > Norm der Funktion\n";
  } catch (DROPSErrCL d) {
        d.handle();
    }

    return 0;
}

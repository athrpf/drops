#include "geom/multigrid.h"
#include "out/output.h"
#include "geom/builder.h"
#include "stokes/stokes.h"
#include "num/stokessolver.h"
#include "out/output.h"
#include <fstream>

// q*u - nu*laplace u + Dp = f
//                  -div u = 0
class DrivenCavityCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static DROPS::SVectorCL<3> f(const DROPS::Point3DCL&, double)
        { DROPS::SVectorCL<3> ret(0.0); return ret; }
    const double nu;
    
    DrivenCavityCL() : nu(1.0) {}
};

DROPS::SVectorCL<3> Null( const DROPS::Point3DCL&, double)   { return DROPS::SVectorCL<3>(0.); }
DROPS::SVectorCL<3> Stroem( const DROPS::Point3DCL&, double) { DROPS::SVectorCL<3> ret(0.); ret[0]= 1.; return ret; }

void MarkLower( DROPS::MultiGridCL& mg, double tresh)
{
    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin()),
             ItEnd(mg.GetTriangTetraEnd()); It!=ItEnd; ++It)
    {
        if (GetBaryCenter(*It)[2]<=tresh )
            It->SetRegRefMark();
    }
}


namespace DROPS // for Strategy
{

template<class Coeff>
void Strategy(StokesP2P1CL<Coeff>& Stokes, double inner_iter_tol, Uint maxStep)
// flow control
{
    MultiGridCL& MG= Stokes.GetMG();

    IdxDescCL  loc_vidx, loc_pidx;
    IdxDescCL* vidx1= &Stokes.vel_idx;
    IdxDescCL* pidx1= &Stokes.pr_idx;
    IdxDescCL* vidx2= &loc_vidx;
    IdxDescCL* pidx2= &loc_pidx;
//    IdxDescCL* err_idx= &_err_idx;
    VecDescCL     loc_p;
    VelVecDescCL  loc_v;
    VelVecDescCL* v1= &Stokes.v;
    VelVecDescCL* v2= &loc_v;
    VecDescCL*    p1= &Stokes.p;
    VecDescCL*    p2= &loc_p;
    VelVecDescCL* b= &Stokes.b;
    VelVecDescCL* c= &Stokes.c;
//    VecDescCL* err= &_err;
    MatDescCL* A= &Stokes.A;
    MatDescCL* B= &Stokes.B;
    Uint step= 0;
//    bool new_marks;
//    double akt_glob_err;

    vidx1->Set( 3, 3, 0, 0);
    vidx2->Set( 3, 3, 0, 0);
    pidx1->Set( 1, 0, 0, 0);
    pidx2->Set( 1, 0, 0, 0);
    
    TimerCL time;
//    err_idx->Set( 0, 0, 0, 1);
    do
    {
//        akt_glob_err= glob_err;
        MarkLower(MG,0.25); 
//        MarkAll(MG); 
        MG.Refine();
        Stokes.CreateNumberingVel(MG.GetLastLevel(), vidx1);    
        Stokes.CreateNumberingPr(MG.GetLastLevel(), pidx1);    
        std::cerr << "altes und neues TriangLevel: " << vidx2->TriangLevel << ", "
                  << vidx1->TriangLevel << std::endl;
        MG.SizeInfo(std::cerr);
        b->SetIdx(vidx1);
        c->SetIdx(pidx1);
        p1->SetIdx(pidx1);
        v1->SetIdx(vidx1);
        std::cerr << "Anzahl der Druck-Unbekannten: " << p2->Data.size() << ", "
                  << p1->Data.size() << std::endl;
        std::cerr << "Anzahl der Geschwindigkeitsunbekannten: " << v2->Data.size() << ", "
                  << v1->Data.size() << std::endl;
        if (p2->RowIdx)
        {
//            const StokesBndDataCL& BndData= Stokes.GetBndData();
//            P1EvalCL<double, const StokesPrBndDataCL, const VecDescCL>  pr2(p2, &BndData.Pr, &MG);
//            P1EvalCL<double, const StokesPrBndDataCL, VecDescCL>        pr1(p1, &BndData.Pr, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, const VelVecDescCL> vel2(v2, &BndData.Vel, &MG);
//            P2EvalCL<SVectorCL<3>, const StokesVelBndDataCL, VelVecDescCL>       vel1(v1, &BndData.Vel, &MG);
//            Interpolate(pr1, pr2);
//            Interpolate(vel1, vel2);
//            CheckSolution(v1,p1,&LsgVel,&LsgPr);
            v2->Reset();
            p2->Reset();
        }
        A->Reset();
        B->Reset();
        A->SetIdx(vidx1, vidx1);
        B->SetIdx(pidx1, vidx1);
        time.Reset();
        time.Start();
        Stokes.SetupSystem(A, b, B, c);
        time.Stop();
        std::cerr << time.GetTime() << " seconds for setting up all systems!" << std::endl;
        time.Reset();
        time.Start();
        A->Data * v1->Data;
        time.Stop();
        std::cerr << " A*x took " << time.GetTime() << " seconds!" << std::endl;
        time.Reset();
        time.Start();
        transp_mul( A->Data, v1->Data);
        time.Stop();
        std::cerr << "AT*x took " << time.GetTime() << " seconds!" << std::endl;
/*        
        { // write system in files for MatLab
            std::ofstream Adat("Amat.dat"), Bdat("Bmat.dat"), bdat("fvec.dat"), cdat("gvec.dat");
            Adat << A->Data;   Bdat << B->Data;    bdat << b->Data;    cdat << c->Data;
        }
*/ //        Stokes.GetDiscError(&LsgVel, &LsgPr);
//std::cout << A->Data << std::endl << b->Data << std::endl;
/*        double half= M_PI/8;
        MultiGridCL::TriangVertexIteratorCL vert= MG.GetTriangVertexBegin(A->GetRowLevel());
        while (vert->GetCoord()[0]!=half || vert->GetCoord()[1]!=half || vert->GetCoord()[2]!=half) ++vert;
        IdxT unk= vert->Unknowns(A->RowIdx->Idx);
        std::cerr << vert->GetCoord() << " has index " << unk << std::endl;
        std::cerr << "A(i,i) = " << A->Data(unk,unk) <<std::endl;    
        std::cerr << "B(i,j) = " << B->Data(vert->Unknowns(B->RowIdx->Idx),unk) << std::endl;
*/        Uint meth;
        std::cerr << "\nwhich method? 0=Uzawa, 1=Schur > "; std::cin >> meth;
        time.Reset();

        MatDescCL M;
        M.SetIdx( pidx1, pidx1);
        Stokes.SetupPrMass( &M);

        double outer_tol;
        std::cerr << "tol = "; std::cin >> outer_tol;

        if (meth)
        {
//            PSchur_PCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
            PSchur_GSPCG_CL schurSolver( M.Data, 200, outer_tol, 200, inner_iter_tol);
            time.Start();
            schurSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
        }
        else // Uzawa
        {
            double tau;
            Uint inner_iter;
            std::cerr << "tau = "; std::cin >> tau;
            std::cerr << "#PCG steps = "; std::cin >> inner_iter;
            Uzawa_PCG_CL uzawaSolver( M.Data, 5000, outer_tol, inner_iter, inner_iter_tol, tau);
            time.Start();
            uzawaSolver.Solve( A->Data, B->Data, v1->Data, p1->Data, b->Data, c->Data);
            time.Stop();
            std::cerr << "Iterationen: " << uzawaSolver.GetIter()
                      << "\tNorm des Res.: " << uzawaSolver.GetResid() << std::endl;
        }
        std::cerr << "Das Verfahren brauchte "<<time.GetTime()<<" Sekunden.\n";
//        Stokes.CheckSolution(v1, p1, &LsgVel, &LsgPr);
        A->Reset();
        B->Reset();
        b->Reset();
        c->Reset();
//        std::cerr << "Loesung Druck: " << p1->Data << std::endl;
//        CreateNumbering(new_idx->TriangLevel, err_idx);
//        err->SetIdx(err_idx);
//        NumMarked= EstimateError(new_x, 0.2, &Estimator);
//TODO: Fehler schaetzen
//        new_marks= EstimateError(new_x, 1.5, akt_glob_err, &ResidualErrEstimator);
//        err->Reset();
//        DeleteNumbering(err_idx);
        std::swap(v2, v1);
        std::swap(p2, p1);
        std::swap(vidx2, vidx1);
        std::swap(pidx2, pidx1);
        std::cerr << std::endl;
    }
    while (++step<maxStep);
    // we want the solution to be in Stokes.v, Stokes.pr
    if (v2 == &loc_v)
    {
        Stokes.vel_idx.swap( loc_vidx);
        Stokes.pr_idx.swap( loc_pidx);
        Stokes.v.SetIdx( &Stokes.vel_idx);
        Stokes.p.SetIdx( &Stokes.pr_idx);
        
        Stokes.v.Data= loc_v.Data;
        Stokes.p.Data= loc_p.Data;
    }
}

} // end of namespace DROPS


void MarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRegRefMark();
    }
}


void UnMarkDrop (DROPS::MultiGridCL& mg, DROPS::Uint maxLevel)
{
    DROPS::Point3DCL Mitte; Mitte[0]=0.5; Mitte[1]=0.5; Mitte[2]=0.5;

    for (DROPS::MultiGridCL::TriangTetraIteratorCL It(mg.GetTriangTetraBegin(maxLevel)),
             ItEnd(mg.GetTriangTetraEnd(maxLevel)); It!=ItEnd; ++It)
    {
        if ( (GetBaryCenter(*It)-Mitte).norm()<=std::max(0.1,1.5*std::pow(It->GetVolume(),1.0/3.0)) )
            It->SetRemoveMark();
    }
}

int main (int argc, char** argv)
{
  try
  {
    if (argc!=3)
    {
        std::cerr << "You have to specify two parameters:\n\tdrivcav <inner_iter_tol> <num_refinement>" << std::endl;
        return 1;
    }

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= e2[1]= e3[2]= 1.;

    typedef DROPS::StokesP2P1CL<DrivenCavityCL> 
            StokesOnBrickCL;
    typedef StokesOnBrickCL MyStokesCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
//    DROPS::BBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2, 2);
//    DROPS::LBuilderCL brick(null, e1, e2, e3, 4, 4, 4, 2, 2);
    const bool IsNeumann[6]= 
        {false, false, false, false, false, false};
    const DROPS::StokesBndDataCL::VelBndDataCL::bnd_val_fun bnd_fun[6]= 
        { &Null, &Null, &Null, &Null, &Null, &Stroem };
        
    StokesOnBrickCL prob(brick, DrivenCavityCL(), DROPS::StokesBndDataCL(6, IsNeumann, bnd_fun));
    DROPS::MultiGridCL& mg = prob.GetMG();
    DROPS::RBColorMapperCL colormap;
//    MarkAll(mg);


//    mg.Refine();
//   MarkAll(mg);
//    mg.Refine();
/*
    for (DROPS::Uint i=0; i<4; ++i)
    {
        MarkDrop(mg, std::min(mg.GetLastLevel(),10u));
//        MarkAll(mg);
//        MarkSome(mg);
        std::cerr << DROPS::SanityMGOutCL(mg);
//        std::cerr << DROPS::DumpMGCL(mg);
//        std::cerr << i << std::endl;
        mg.Refine();
    }
*/
//    MarkAll(mg);
    double inner_iter_tol= atof(argv[1]);
    int num_ref= atoi(argv[2]);
    std::cerr << "inner iter tol: " << inner_iter_tol << std::endl;
    std::cerr << "refinements: " << num_ref << std::endl;
    Strategy(prob, inner_iter_tol, num_ref);
    std::cerr << DROPS::SanityMGOutCL(mg) << std::endl;
    std::ofstream fil("ttt.off");
    double min= prob.p.Data.min(),
           max= prob.p.Data.max();
    std::cerr << "pressure min/max: "<<min<<", "<<max<<std::endl;
    fil << DROPS::GeomSolOutCL<MyStokesCL::const_DiscPrSolCL>(mg, prob.GetPrSolution(), &colormap, -1, false, 0.0, -10, 10) << std::endl;

    DROPS::IdxDescCL tecIdx;
    tecIdx.Set( 1, 0, 0, 0);
    prob.CreateNumberingPr( mg.GetLastLevel(), &tecIdx);    
std::cerr <<"TecIdx = "<< tecIdx.GetIdx()<<std::endl;    
    std::ofstream v2d("data2D.dat");
    DROPS::TecPlot2DSolOutCL< MyStokesCL::const_DiscVelSolCL, MyStokesCL::const_DiscPrSolCL>
        tecplot2d( mg, prob.GetVelSolution(), prob.GetPrSolution(), tecIdx, -1, 1, 0.5); // cutplane is y=0.5
    v2d << tecplot2d;
    v2d.close();

    std::ofstream v3d("data3D.dat");
//    DROPS::TecPlotSolOutCL< MyStokesCL::const_DiscVelSolCL, MyStokesCL::const_DiscPrSolCL>
//        tecplot( mg, prob.GetVelSolution(), prob.GetPrSolution() );
//    v3d << tecplot;
    DROPS::MapleSolOutCL<MyStokesCL::const_DiscVelSolCL, MyStokesCL::const_DiscPrSolCL>
        mapleplot( mg, prob.GetVelSolution(), prob.GetPrSolution(), -1, DROPS::PlaneCL(DROPS::std_basis<3>(2), .5));
    v3d << mapleplot;
    v3d.close();

/*
    int nv= 2*cbrt( mg.GetVertices().size()) - 1; // cbrt = 3. Wurzel;
    int n= nv-1;  // = Anzahl der Bloecke in einer Dimension
    
    typedef void* voidpointerT;
    void **data= new voidpointerT[nv*nv*nv];
    DROPS::Point3DCL coord;
    for (DROPS::MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(), end= mg.GetTriangVertexEnd(); it!=end; ++it)
    {
        coord= it->GetCoord();
        int idx= rint( coord[0]*n)*nv*nv + rint( coord[1]*n)*nv + rint( coord[2]*n);
        data[idx]= static_cast<voidpointerT>(&*it);
    }
    for (DROPS::MultiGridCL::TriangEdgeIteratorCL it= mg.GetTriangEdgeBegin(), end= mg.GetTriangEdgeEnd(); it!=end; ++it)
    {
        coord= DROPS::GetBaryCenter( *it);
        int idx= rint( coord[0]*n)*nv*nv + rint( coord[1]*n)*nv + rint( coord[2]*n);
        data[idx]= static_cast<voidpointerT>(&*it);
    }
    MyStokesCL::DiscVelSolCL vel= prob.GetVelSolution();
    MyStokesCL::DiscPrSolCL pr= prob.GetPrSolution();
    std::ofstream v2d("data2D.dat");
    v2d << " VARIABLES = " << '"' << 'X' << '"' << " , "
                           << '"' << 'Z' << '"' << " , "
                           << '"' << 'U' << '"' << " , "
                           << '"' << 'W' << '"' << " , "
                           << '"' << 'P' << '"' << std::endl;
    v2d << "ZONE I = "<<nv<<", J = "<<nv<<", F = POINT"<<std::endl;
    const int ny= nv/2;
    for (int nz=0; nz<nv; ++nz)
        for (int nx=0; nx<nv; ++nx)
        {
            if (nx%2 || nz%2) // nx or nz is odd
            {
                DROPS::EdgeCL* ep= reinterpret_cast<DROPS::EdgeCL*>(data[nx*nv*nv + ny*nv + nz]);
                coord= DROPS::GetBaryCenter( *ep);
                v2d << coord[0] << "  " << coord[2] << "  "
                    << vel.val( *ep)[0] << "  " << vel.val( *ep)[2] << "  " << pr.val( *ep, 0.5) <<'\n';
            }
            else
            {
                DROPS::VertexCL* vp= reinterpret_cast<DROPS::VertexCL*>(data[nx*nv*nv + ny*nv + nz]);
                coord= vp->GetCoord();
                v2d << coord[0] << "  " << coord[2] << "  "
                    << vel.val( *vp)[0] << "  " << vel.val( *vp)[2] << "  " << pr.val( *vp) <<'\n';
            }
        }
        
    v2d.close();
    n= (nv+1)/2;
    std::ofstream v3d("data3D.dat");
    v3d << " VARIABLES = " << '"' << 'X' << '"' << " , "
                           << '"' << 'Y' << '"' << " , "
                           << '"' << 'Z' << '"' << " , "
                           << '"' << 'U' << '"' << " , "
                           << '"' << 'V' << '"' << " , "
                           << '"' << 'W' << '"' << " , "
                           << '"' << 'P' << '"' << std::endl;
    v3d << "ZONE I = "<<n<<", J = "<<n<<", K = "<<n<<", F = POINT"<<std::endl;
    for (DROPS::MultiGridCL::TriangVertexIteratorCL it= mg.GetTriangVertexBegin(), end= mg.GetTriangVertexEnd(); it!=end; ++it)
        {
            v3d << it->GetCoord()[0] << "  " << it->GetCoord()[1] << "  " << it->GetCoord()[2] << "  "
                << vel.val( *it)[0] << "  " << vel.val( *it)[1] << "  " << vel.val( *it)[2] << "  " << pr.val( *it)<< '\n';
        }
    v3d.close();
*/
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
}

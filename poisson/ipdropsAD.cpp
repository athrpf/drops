//**************************************************************************
// File:    ipdropsAD.cpp                                                  *
// Content: ipdrops, test case for AD                                      *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - May, 02 2006                                           *
//**************************************************************************

#include "poisson/instatpoisson.h"
#include "num/solver.h"
#include "poisson/integrTime.h"
#include "out/output.h"
#include <fstream>


// Das folgende Problem wird untersucht
// du/dt - nu * laplace(u) + q*u = f

class ParamsCL // Parameter des Problems
{
  public:
    double nu, dt, xl, yl, zl, M, cgtol, theta, q_h;
    int cgiter, nx, ny, nz, nt;

    ParamsCL()
    {
        M = 0; // Rueckgabe der Temperatur fuer x=M
        xl= 0.025; // Abmessung des Gebietes
        yl= 19.5;
        zl= 39;
        nx = 2; // Anzahl der Intervalle
        ny = 4;
        nz = 20;
        dt = 2e-3;
        nt = 100;
        nu= 6.303096; // Temperaturleitfaehigkeit
        q_h= 278.652e-3; // Heizwaermestrom
        theta = 0.5; // Crank-Nicholson
        cgtol = 1e-6;
        cgiter = 200;
    }
} C;


// Setzen der Koeffizienten
class PoissonCoeffCL
{
  public:
    // static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double alpha(const DROPS::Point3DCL&, double)
      { return 1; }
    static double f(const DROPS::Point3DCL&, double)
    {
        return 0.0;
    }
    static DROPS::Point3DCL Vel(const DROPS::Point3DCL&, double)
      { return DROPS::Point3DCL(0.); } // no convection
};

double SinusWave( const DROPS::Point3DCL& p, double t)
{
    return C.q_h + 200e-3*sin((p[2]/C.zl+t/0.2)*4*M_PI);
}


// Die Klasse MatConnect bildet das Interface zwischen
// C++ und Matlab

class MatConnect
{
  private:
    typedef std::pair<double, double> d_pair;
    typedef std::pair<double, d_pair> cmp_key;
    typedef std::map<cmp_key, double*> node_map;
    typedef node_map::const_iterator ci;
    typedef DROPS::InstatPoissonBndDataCL::bnd_val_fun BndFuncT;

    node_map _NodeMap;
    static const DROPS::VectorCL* _BndData[6];
    const DROPS::VectorCL* _T0;
    DROPS::VectorCL _weight;
    static double _DeltaT, _XLen, _YLen, _ZLen, _SpIncrX, _SpIncrY, _SpIncrZ;
    double _CutPos;
    static int _nt, _MeshRefX, _MeshRefY, _MeshRefZ, _FacePtsYZ;
    int _Count;

  public:
    MatConnect(double DeltaT, int time_steps,
      double xl, double yl, double zl,
      int mrx, int mry, int mrz,
      double CutPos, int FacePtsYZ,
      const DROPS::VectorCL& T0,
      const DROPS::VectorCL& S1, const DROPS::VectorCL& S2, const DROPS::VectorCL& S3,
      const DROPS::VectorCL& S4, const DROPS::VectorCL& S5, const DROPS::VectorCL& S6)
    {
      _NodeMap.clear();
      _DeltaT= DeltaT; _nt= time_steps;
      _XLen= xl; _YLen= yl; _ZLen= zl;
      _MeshRefX= mrx; _MeshRefY= mry; _MeshRefZ= mrz;
      _SpIncrX= xl/mrx; _SpIncrY= yl/mry; _SpIncrZ= zl/mrz;
      _Count= 0; _CutPos= CutPos; _FacePtsYZ= FacePtsYZ;
      _T0= &T0;
      _BndData[0]= &S1; _BndData[1]= &S2; _BndData[2]= &S3;
      _BndData[3]= &S4; _BndData[4]= &S5; _BndData[5]= &S6;

      // init weights for quadrature on cartesian grid
      const int N= FacePtsYZ*(time_steps+1);
      _weight.resize(N);
      for (int i=0; i<N; ++i) _weight[i]= 1;
      // spatial/temporal boundaries get factor 1/2
      for (int i=0; i<FacePtsYZ; ++i)
      {
        _weight[i]*= 0.5;               // t=0
        _weight[N-FacePtsYZ+i]*= 0.5;   // t=nt*dt
      }
      for (int it=0; it<=_nt; ++it)
        for (int iz=0; iz<=_MeshRefZ; ++iz)
        {
          _weight[it*FacePtsYZ+iz*(_MeshRefY+1)]*= 0.5;                 // y=0
          _weight[it*FacePtsYZ+iz*(_MeshRefY+1)+_MeshRefY]*= 0.5;       // y=yl
        }
      for (int it=0; it<=_nt; ++it)
        for (int iy=0; iy<=_MeshRefY; ++iy)
        {
          _weight[it*FacePtsYZ+                         +iy]*= 0.5;     // z=0
          _weight[it*FacePtsYZ+(_MeshRefZ)*(_MeshRefY+1)+iy]*= 0.5;     // z=zl
        }
      std::cerr << "Quadrature weights initialized!\n";
    }

    template<int num> static int getLexNum(const DROPS::Point3DCL& p)
    {
      switch(num)
      {
        case 0: // y-z plane
        {
          int nY= static_cast<int>(p[1]/_SpIncrY+0.5),  nZ= static_cast<int>(p[2]/_SpIncrZ+0.5);
          return nZ*(_MeshRefY+1) + nY;
        }
        case 1: // x-z plane
        {
          int nX= static_cast<int>(p[0]/_SpIncrX+0.5),  nZ= static_cast<int>(p[2]/_SpIncrZ+0.5);
          return nZ*(_MeshRefX+1) + nX;
        }
        case 2: default: // x-y plane
        {
          int nX= static_cast<int>(p[0]/_SpIncrX+0.5),  nY= static_cast<int>(p[1]/_SpIncrY+0.5);
          return nY*(_MeshRefX+1) + nX;
        }
      }
    }

    template<int num> static double getBndVal(const DROPS::Point3DCL& p, double t)
    {
      int count= 0;
      const int nT= static_cast<int>(t/_DeltaT+0.5);
      switch(num)
      {
        case 0: case 1: // y-z plane
          count= nT*(_MeshRefY+1)*(_MeshRefZ+1) + getLexNum<0>(p);
          break;
        case 2: case 3: // x-z plane
          count= nT*(_MeshRefX+1)*(_MeshRefZ+1) + getLexNum<1>(p);
          break;
        case 4: case 5: default: // x-y plane
          count= nT*(_MeshRefX+1)*(_MeshRefY+1) + getLexNum<2>(p);
      }

      return (*_BndData[num])[count];
    }

    template<int num> static void InitVec( DROPS::VectorCL& v, BndFuncT f)
    {
        double t= 0;
        DROPS::Point3DCL p;
        int n0, n1, c0, c1;
        double d0, d1;
        switch (num)
        {
          case 0: case 1: // y-z plane
            c0= 1;            c1= 2;
            n0= _MeshRefY;    n1= _MeshRefZ;
            d0= _SpIncrY;     d1= _SpIncrZ;
            break;
          case 2: case 3: // x-z plane
            c0= 0;            c1= 2;
            n0= _MeshRefX;    n1= _MeshRefZ;
            d0= _SpIncrX;     d1= _SpIncrZ;
            break;
          case 4: case 5: default: // x-y plane
            c0= 0;            c1= 1;
            n0= _MeshRefX;    n1= _MeshRefY;
            d0= _SpIncrX;     d1= _SpIncrY;
        }
        const int FacePts= (n0+1)*(n1+1);

        for (int it= 0; it<=_nt; ++it)
        {
            for (int i1= 0; i1<=n1; ++i1)
            {
                p[c1]= i1*d1;
                for (int i0= 0; i0<=n0; ++i0)
                {
                    p[c0]= i0*d0;
                    const int idx= it*FacePts + getLexNum<num>(p);
                    v[idx]= f(p,t);
                }
            }
            t+= _DeltaT;
        }
    }

    double L2ScalarProd( const DROPS::VectorCL& u, const DROPS::VectorCL& v)
    {
        double sum= 0;
        for (int i=0, n= _FacePtsYZ*(_nt+1); i<n; ++i)
        {
            sum+= _weight[i]*u[i]*v[i];
        }
        return sum*_SpIncrX*_SpIncrY*_SpIncrZ*_DeltaT;
    }

    static void setBndData( int num, const DROPS::VectorCL& S) { _BndData[num]= &S; }
    void setInitialData( const DROPS::VectorCL& S) { _T0= &S; }

    void setNodeMap(DROPS::VecDescCL& x, DROPS::MultiGridCL& MG)
    {
      DROPS::Point3DCL pt;
      DROPS::Uint lvl= x.GetLevel();
      DROPS::Uint indx= x.RowIdx->GetIdx();

      d_pair help;
      cmp_key key;
      for (DROPS::MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(lvl),
        send=MG.GetTriangVertexEnd(lvl); sit != send; ++sit)
      {
        DROPS::IdxT i= sit->Unknowns(indx);
        pt= sit->GetCoord();

        help= std::make_pair(pt[2], pt[1]);
        key= std::make_pair(pt[0], help);
        _NodeMap[key]= &(x.Data[i]);
      }
    }

    void getInitialData()
    {
      double xc, yc, zc;
      for (ci p= _NodeMap.begin(); p!= _NodeMap.end(); p++)
      {
        xc= p->first.first;
        yc= p->first.second.second;
        zc= p->first.second.first;

        const int a= static_cast<int>((_MeshRefX*xc/_XLen)+0.5);
        const int b= static_cast<int>((_MeshRefY*yc/_YLen)+0.5);
        const int c= static_cast<int>((_MeshRefZ*zc/_ZLen)+0.5);

        const int count= a*_FacePtsYZ + c*(_MeshRefY+1) + b;
        //printf("xc, yc, zc, uh: %g %g %g %g\n", xc, yc, zc, (*_T0)[count]);
        //printf("a, b, c: %d %d %d\n", a, b, c);
        //printf("count: %d\n", count);

        *(p->second)= (*_T0)[count];
      }
    }

    void setOutputData(DROPS::VectorCL& sol2D, bool firstStep= false)
    {
      if (firstStep) _Count= 0;
      ci p= _NodeMap.begin();
      double xc= p->first.first;
      while ((xc<_CutPos)&&(p!=_NodeMap.end()))
      {
        p++;
        xc= p->first.first;
      }
      //for (int i=0; i<_CutPos; i++)
      //  for (int j=0; j<_FacePtsYZ; j++)
      //    if (p!=_NodeMap.end())
      //      p++;
      for (int k=0; k<_FacePtsYZ; k++)
      {
        sol2D[_Count]= *(p->second);
        //printf("xc yc zc sol2D: %g %g %g %g\n", p->first.first, p->first.second.second,
        //  p->first.second.first, *(p->second));
        p++;
        _Count++;
      }
    }

    void printData() const
    {
      printf("DeltaT: %g\n", _DeltaT);
      printf("XLen: %g\n", _XLen);
      printf("YLen: %g\n", _YLen);
      printf("ZLen: %g\n", _ZLen);
      printf("SpIncrX: %g\n", _SpIncrX);
      printf("SpIncrY: %g\n", _SpIncrY);
      printf("SpIncrZ: %g\n", _SpIncrZ);
      printf("MeshRefX: %d\n", _MeshRefX);
      printf("MeshRefY: %d\n", _MeshRefY);
      printf("MeshRefZ: %d\n", _MeshRefZ);
      printf("Count: %d\n", _Count);
      printf("CutPos: %g\n", _CutPos);
      printf("FacePtsYZ: %d\n", _FacePtsYZ);

      for (ci p=_NodeMap.begin(); p!=_NodeMap.end(); p++)
        printf("xc yc zc uh: %g %g %g %g\n", p->first.first, p->first.second.second,
          p->first.second.first, *(p->second));
    }
};

const DROPS::VectorCL* MatConnect::_BndData[6];
double MatConnect::_DeltaT= .0;
double MatConnect::_XLen= .0;
double MatConnect::_YLen= .0;
double MatConnect::_ZLen= .0;
double MatConnect::_SpIncrX= .0;
double MatConnect::_SpIncrY= .0;
double MatConnect::_SpIncrZ= .0;
int MatConnect::_nt= 0;
int MatConnect::_MeshRefX= 0;
int MatConnect::_MeshRefY= 0;
int MatConnect::_MeshRefZ= 0;
int MatConnect::_FacePtsYZ= 0;



namespace DROPS // for Strategy
{

double getIsolatedBndVal(const Point3DCL&, double) { return 0.0; }

template<class Coeff>
void Strategy(InstatPoissonP1CL<Coeff>& Poisson, DROPS::VectorCL& sol2D, MatConnect& MatCon)
{
  typedef InstatPoissonP1CL<Coeff> MyPoissonCL;

  MLIdxDescCL& idx= Poisson.idx;
  VecDescCL& x= Poisson.x;
  VecDescCL& b= Poisson.b;
  MLMatDescCL& A= Poisson.A;
  MLMatDescCL& M= Poisson.M;

  VecDescCL cplA;
  VecDescCL cplM;

  idx.Set(1, 0, 0, 0);

  MultiGridCL& MG= Poisson.GetMG();

  // erzeuge Nummerierung zu diesem Index
  Poisson.CreateNumbering(MG.GetLastLevel(), &idx);

  // Vektoren mit Index idx
  b.SetIdx(&idx);
  x.SetIdx(&idx);
  cplA.SetIdx(&idx);
  cplM.SetIdx(&idx);

  printf("Anzahl der Unbekannten: %d\n", x.Data.size());
//  printf("Theta: %g\n", C.theta);
//  printf("Toleranz CG: %g\n", C.cgtol);
//  printf("max. Anzahl CG-Iterationen: %d\n", C.cgiter);

  // Steifigkeitsmatrix mit Index idx (Zeilen und Spalten)
  A.SetIdx(&idx, &idx);
  // Massematrix mit Index idx (Zeilen und Spalten)
  M.SetIdx(&idx, &idx);

  // stationaerer Anteil
  Poisson.SetupInstatSystem(A, M, Poisson.t);

  // instationaere rechte Seite
  Poisson.SetupInstatRhs(cplA, cplM, Poisson.t, b, Poisson.t);

  // PCG-Verfahren mit SSOR-Vorkonditionierer
  SSORPcCL pc(1.0);
  PCG_SsorCL pcg_solver(pc, C.cgiter, C.cgtol);

  // Zeitdiskretisierung mit one-step-theta-scheme
  // theta=1 -> impl. Euler
  // theta=0.5 -> Crank-Nicholson
  InstatPoissonThetaSchemeCL<InstatPoissonP1CL<Coeff>, PCG_SsorCL>
    ThetaScheme(Poisson, pcg_solver, C.theta);
  ThetaScheme.SetTimeStep(C.dt, C.nu);
  MatCon.setNodeMap(x, MG);
  MatCon.getInitialData();
  MatCon.setOutputData(sol2D, true);

  int MaxIterCG= 0;
  double MaxResCG= 0.;

  for (int step=1;step<=C.nt;step++)
  {
    ThetaScheme.DoStep(x);
    if (MaxIterCG<=pcg_solver.GetIter())
      MaxIterCG= pcg_solver.GetIter();
    if (MaxResCG<=pcg_solver.GetResid())
      MaxResCG= pcg_solver.GetResid();

    //printf("t= %g\n", Poisson.t);
    //printf("Iterationen: %d", pcg_solver.GetIter());
    //printf("    Norm des Residuums: %g\n", pcg_solver.GetResid());

    MatCon.setOutputData(sol2D);
  }

  printf("\nAnzahl CG-Iterationen (max. pro Zeitschritt): %d\n", MaxIterCG);
  printf("Norm des max. Residuums CG: %g\n", MaxResCG);

  //MatCon.printData();

  A.Reset();
  b.Reset();
}

} // end of namespace DROPS


static
void ipdrops(DROPS::VectorCL& sol2D, MatConnect& MatCon)
{
  try
  {
    //DROPS::TimerCL time;
    //time.Reset();

    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= C.xl;
    e2[1]= C.yl;
    e3[2]= C.zl;

    typedef DROPS::InstatPoissonP1CL<PoissonCoeffCL>
      InstatPoissonOnBrickCL;
    typedef InstatPoissonOnBrickCL MyPoissonCL;

    DROPS::BrickBuilderCL brick(null, e1, e2, e3, C.nx, C.ny, C.nz);

    printf("\nDelta t = %g", C.dt);
    printf("\nAnzahl der Zeitschritte = %d\n", C.nt);

    // Randdaten
    const bool isneumann[6]= { true, true, true, true, true, true };
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
      { &MatConnect::getBndVal<0>, &MatConnect::getBndVal<1>,
        &DROPS::getIsolatedBndVal, &DROPS::getIsolatedBndVal,
        &DROPS::getIsolatedBndVal, &DROPS::getIsolatedBndVal };

    DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);

    // mg.SizeInfo();
    DROPS::Strategy(prob, sol2D, MatCon);

    //time.Stop();
    //printf("Zeit fuer das Loesen des direkten Problems: %g sek\n", time.GetTime());

    return;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}

double Guetefunktional( DROPS::VectorCL& q_Film)
{
    /* Set the vectors. */
    const int DiscPtsXY= (C.nx+1)*(C.ny+1)*(C.nt+1);
    const int DiscPtsXZ= (C.nx+1)*(C.nz+1)*(C.nt+1);
    const int DiscPtsYZ= (C.ny+1)*(C.nz+1)*(C.nt+1);

    DROPS::VectorCL T0( 36.8097, (C.nx+1)*(C.ny+1)*(C.nz+1)); // mittlere Fluidtemperatur im Gesamtgebiet fuer Zeitschritt 0
    DROPS::VectorCL q_Film_exakt(DiscPtsYZ), T_mess(DiscPtsYZ), sol2D(DiscPtsYZ),
        q_heat( -C.q_h, DiscPtsYZ); // Heizwaermestrom


    DROPS::VectorCL VecXY(DiscPtsXY); // Nullvektor
    DROPS::VectorCL VecXZ(DiscPtsXZ);

    MatConnect MatCon(C.dt, C.nt, C.xl, C.yl, C.zl, C.nx, C.ny, C.nz, C.M, (C.ny+1)*(C.nz+1),
      T0, q_heat, q_Film_exakt, VecXZ, VecXZ, VecXY, VecXY);

    MatCon.InitVec<1>( q_Film_exakt, &SinusWave);

    /* Call the solver. */
    // erzeuge Messdaten zu exaktem q_Film
    MatCon.setBndData( 1, q_Film_exakt);
    ipdrops( T_mess, MatCon);

    // Temperatur zu uebergebenem q_Film
    MatCon.setBndData( 1, q_Film);
    ipdrops( sol2D, MatCon);

    DROPS::VectorCL diff( T_mess - sol2D);
    const double J= MatCon.L2ScalarProd( diff, diff);
    return J;
}


int main()
{
    DROPS::VectorCL q( C.q_h, (C.ny+1)*(C.nz+1)*(C.nt+1));
    const double J= Guetefunktional( q);
    std::cerr << "J = " << J << std::endl;
    return 0;
}

//**************************************************************************
// File:    ipdrops.cpp                                                    *
// Content: program with interface for matlab                              *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - Nov, 19 2002                                           *
//**************************************************************************

#include "mex.h"

#include "../instatpoisson.h"
#include "../../num/solver.h"
#include "../integrTime.h"


// Das folgende Problem wird untersucht
// du/dt - nu * laplace(u) + q*u = f


extern void _main();


// Setzen der Koeffizienten
class PoissonCoeffCL
{	
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double f(const DROPS::Point3DCL& p, double t) { return 0.0; }
};


// Die Klasse MatConnect bildet das Interface zwischen
// C++ und Matlab

class MatConnect
{
  private:
    static double* _BndData[6];
    static double _DeltaT, _XLen, _YLen, _ZLen, _SpIncrX, _SpIncrY, _SpIncrZ;
    static int _MeshRefX, _MeshRefY, _MeshRefZ;
	
  public:
    MatConnect(double DeltaT, 
      double xl, double yl, double zl,
      int mrx, int mry, int mrz,  
      double* S1, double* S2, double* S3,
      double* S4, double* S5, double* S6)
    {
      _DeltaT= DeltaT;
      _XLen= xl; _YLen= yl; _ZLen= zl;
      _MeshRefX= mrx; _MeshRefY= mry; _MeshRefZ= mrz;
      _SpIncrX= xl/mrx; _SpIncrY= yl/mry; _SpIncrZ= zl/mrz;
      _BndData[0]= S1; _BndData[1]= S2; _BndData[2]= S3;
      _BndData[3]= S4; _BndData[4]= S5; _BndData[5]= S6;
    }
    
    template<int num> static int getLexNum(const DROPS::Point2DCL& p)
    {
      int count= 0;
      if (num==0)  // Punkt liegt auf der Seite S1 oder S2
      {
        for (double zcoord= .0; zcoord<= _ZLen*p[1]-_SpIncrZ/2.0; zcoord+= _SpIncrZ)
          count+= (_MeshRefY+1);
        for (double ycoord= .0; ycoord<= _YLen*p[0]+_SpIncrY/2.0; ycoord+= _SpIncrY)
          count++;
      }
      else if (num==1)  // Punkt liegt auf einer der Seiten S3 oder S4
      {
        for (double zcoord= .0; zcoord<= _ZLen*p[1]-_SpIncrZ/2.0; zcoord+= _SpIncrZ)
          count+= (_MeshRefX+1);
        for (double xcoord= .0; xcoord<= _XLen*p[0]+_SpIncrX/2.0; xcoord+= _SpIncrX)
          count++;
      }
      else  // Punkt liegt auf einer der Seiten S5 oder S6
      {
        for (double ycoord= .0; ycoord<= _YLen*p[1]-_SpIncrY/2.0; ycoord+= _SpIncrY)
          count+= (_MeshRefX+1);
        for (double xcoord= .0; xcoord<= _XLen*p[0]+_SpIncrX/2.0; xcoord+= _SpIncrX)
          count++;
      }
    
      return count;
    }
  
    template<int num> static double getBndVal(const DROPS::Point2DCL& p, double t)
    {
      int count= -1;  // das Feld beginnt mit Index 0
    
      // die Matrizen S1..S6 werden spaltenweise uebergeben
      if (num==0 || num==1)  // Seite S1 oder S2
      {
        for (double time= .0; time< t-_DeltaT/2.0; time+= _DeltaT)
          count+= (_MeshRefY+1)*(_MeshRefZ+1);
        count+= getLexNum<0>(p);
      }
      else if (num==2 || num==3) // Seite S3 oder S4
      {
        for (double time= .0; time< t-_DeltaT/2.0; time+= _DeltaT)
          count+= (_MeshRefX+1)*(_MeshRefZ+1);
        count+= getLexNum<1>(p);
      }
      else  // Seite S5 oder S6
      {
        for (double time= .0; time< t-_DeltaT/2.0; time+= _DeltaT)
          count+= (_MeshRefX+1)*(_MeshRefY+1);
        count+= getLexNum<2>(p);
      }
    
      return *(_BndData[num]+count);
    }
};

double* MatConnect::_BndData[6];
double MatConnect::_DeltaT= .0;
double MatConnect::_XLen= .0;
double MatConnect::_YLen= .0;
double MatConnect::_ZLen= .0;
double MatConnect::_SpIncrX= .0;
double MatConnect::_SpIncrY= .0;
double MatConnect::_SpIncrZ= .0;
int MatConnect::_MeshRefX= 0;
int MatConnect::_MeshRefY= 0;
int MatConnect::_MeshRefZ= 0;


namespace DROPS // for Strategy
{

template<class MGB, class Coeff>
void Strategy(InstatPoissonP1CL<MGB, Coeff>& Poisson, double* sol2D, 
  double* T0, int cut_pos, int face_pts, double nu, double dt, 
  int time_steps, double theta, double cgtol, int cgiter)
{
  typedef InstatPoissonP1CL<MGB,Coeff> MyPoissonCL;
  
  IdxDescCL& idx= Poisson.idx;
  VecDescCL& x= Poisson.x;
  VecDescCL& b= Poisson.b;
  MatDescCL& A= Poisson.A;
  MatDescCL& M= Poisson.M;
  
  VecDescCL cplA;
  VecDescCL cplM;
  
  // aktueller Zeitpunkt
  double t= 0;
   
  idx.Set(0, 1, 0, 0, 0);
  
  MultiGridCL& MG= Poisson.GetMG();
  
  // erzeuge Nummerierung zu diesem Index
  Poisson.CreateNumbering(MG.GetLastLevel(), &idx);
  
  // Vektoren mit Index idx
  b.SetIdx(&idx);
  x.SetIdx(&idx);
  cplA.SetIdx(&idx);
  cplM.SetIdx(&idx);
  
  mexPrintf("Anzahl der Unbekannten: %d\n", x.Data.size());
  mexPrintf("Theta: %g\n", theta);
  mexPrintf("Toleranz CG: %g\n", cgtol);
  mexPrintf("max. Anzahl CG-Iterationen: %d\n", cgiter);
  
  // Steifigkeitsmatrix mit Index idx (Zeilen und Spalten)
  A.SetIdx(&idx, &idx);
  // Massematrix mit Index idx (Zeilen und Spalten)
  M.SetIdx(&idx, &idx);
  
  // stationaerer Anteil
  Poisson.SetupInstatSystem(A, M);
  
  // instationaere rechte Seite
  Poisson.SetupInstatRhs(cplA, cplM, t, b, t);
  
  // PCG-Verfahren mit SSOR-Vorkonditionierer
  SSORPcCL pc(1.0);
  PCG_SsorCL pcg_solver(pc, cgiter, cgtol);
  
  // Zeitdiskretisierung mit one-step-theta-scheme
  // theta=1 -> impl. Euler
  // theta=0.5 -> Crank-Nicholson
  InstatPoissonThetaSchemeCL<InstatPoissonP1CL<MGB, Coeff>, PCG_SsorCL>
    ThetaScheme(Poisson, pcg_solver, theta);
  ThetaScheme.SetTimeStep(dt, nu);
  
  
  typedef std::pair<double, double> d_pair;
  typedef std::pair<double, d_pair> cmp_key;
  typedef std::map<cmp_key, double*> node_map;
  typedef node_map::const_iterator ci;
  
  Point3DCL pt;
  Uint lvl= x.RowIdx->TriangLevel;
  Uint indx= x.RowIdx->Idx;
  
  d_pair help;
  cmp_key key;
  node_map nmap;
  int count= 0;
  
  for (MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(lvl),
    send=MG.GetTriangVertexEnd(lvl); sit != send; ++sit)
  {
    IdxT i= sit->Unknowns(indx)[0];
    pt= sit->GetCoord();
    
    help= std::make_pair(pt[2], pt[1]);
    key= std::make_pair(pt[0], help);
    nmap[key]= &(x.Data[i]);
  }
  
  // setze Startwert
  for (ci p= nmap.begin(); p!= nmap.end(); p++)
  {
    *(p->second)= *(T0+count);
    // mexPrintf("%g",*(T0+count));
    count++;
  }


  count= 0;
  for (int step=1;step<=time_steps;step++)
  {
    ThetaScheme.DoStep(x);
    //mexPrintf("t= %g\n", Poisson.t);
    //mexPrintf("Iterationen: %d", pcg_solver.GetIter());
    //mexPrintf("    Norm des Residuums: %g\n", pcg_solver.GetResid());
    //Poisson.CheckSolution(exact_sol, Poisson.t);
    
    // Aufbereitung der Ausgabedaten
    ci p= nmap.begin();
    for (int i=0; i<cut_pos; i++)
      for (int j=0; j<face_pts; j++)
        if (p!=nmap.end())
          p++;
    for (int k=0; k<face_pts; k++)
    {
      //mexPrintf("%g %g %g %g\n", p->first.first, p->first.second.second,
      //  p->first.second.first, *(p->second));
      *(sol2D+count)= *(p->second);
      p++;
      count++;
    }
        
    //for (ci p=nmap.begin(); p!=nmap.end(); p++)
    //  mexPrintf("%g %g %g\n", p->first.first, p->first.second.second,
    //	  p->first.second.first);
  }
  
  
  A.Reset();
  b.Reset();
  
  /*
  // Ausgabe Loesung   
  for (ci p= nmap.begin(); p!= nmap.end(); p++)
  {
    std::cerr << *(p->second) << "\n";
  }
  */
  
}

} // end of namespace DROPS


static
void ipdrops(double* sol2D, double* T0, double* S1, double* S2, 
  double M, double xl, double yl, double zl, double nu, double mrx, 
  double mry, double mrz, double dt, int time_steps, double theta, 
  double cgtol, double cgiter)
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= xl;
    e2[1]= yl;
    e3[2]= zl;
    
    typedef DROPS::InstatPoissonP1CL<DROPS::BrickBuilderCL, PoissonCoeffCL>
      InstatPoissonOnBrickCL;
    typedef InstatPoissonOnBrickCL MyPoissonCL;
    
    int imrx= static_cast<int>(mrx+0.5);
    int imry= static_cast<int>(mry+0.5);
    int imrz= static_cast<int>(mrz+0.5);
    int iM= static_cast<int>(M+0.5);
    int icgiter= static_cast<int>(cgiter+0.5);
    
    /*
    mexPrintf("\nmrx = %d", imrx);
    mexPrintf("\nmry = %d", imry);
    mexPrintf("\nmrz = %d", imrz);  
    */
    
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, imrx, imry, imrz);
    
    mexPrintf("\nRueckgabe der Daten fuer die Flaeche x=%g\n", 
      M*xl/imrx);
    mexPrintf("\nDelta t = %g", dt);
    mexPrintf("\nAnzahl der Zeitschritte = %d\n", time_steps);
    
    // Randdaten 
    int DiscPtsXY= (imrx+1)*(imry+1)*(time_steps+1);
    int DiscPtsXZ= (imrx+1)*(imrz+1)*(time_steps+1);
    DROPS::VectorCL VecXY(DiscPtsXY);
    DROPS::VectorCL VecXZ(DiscPtsXZ);
    double* Sy= &VecXY[0];
    double* Sz= &VecXZ[0];
    
    MatConnect MatCon(dt, xl, yl, zl, imrx, imry, imrz, 
      S1, S2, Sy, Sy, Sz, Sz);  
    
    const bool isneumann[6]= { true, true, true, true, true, true };
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
      { &MatConnect::getBndVal<0>, &MatConnect::getBndVal<1>,
        &MatConnect::getBndVal<2>, &MatConnect::getBndVal<3>,
        &MatConnect::getBndVal<4>, &MatConnect::getBndVal<5> };
    
    DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();
    
    // mg.SizeInfo();
    int FacePtsYZ= (imry+1)*(imrz+1); 
    DROPS::Strategy(prob, sol2D, T0, iM, FacePtsYZ, nu, dt, time_steps,
      theta, cgtol, icgiter);
    
    return;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double nu, dt, xl, yl, zl, M, mrx, mry, mrz, CGTol, CGIter, theta;
  double *T0, *S1, *S2, *sol2D;
  int mrows, ncols;
  
  /* Check for proper number of arguments. */
  if(nrhs!=15)
    mexErrMsgTxt("(T0, S1, S2, M, xl, yl, zl, nu, mrx, mry, mrz, dt, Theta, CGTol, CGIter) as input required.");
  if(nlhs!=1)
    mexErrMsgTxt("Solution on 2D-area as output required.");
  
  /* Check to make sure the first input arguments are double matrices. */
  for (int index=0; index<3; index++)
    if(mxGetPi(prhs[index])!=NULL)
    {
      switch(index) {
      case 0:
        mexErrMsgTxt("Input T0 must be a double matrix.");
      case 1:
        mexErrMsgTxt("Input S1 must be a double matrix.");
      case 2:
        mexErrMsgTxt("Input S2 must be a double matrix.");
      default:
        mexErrMsgTxt("Input error.");
      }
    }
  
  /* Check to make sure the last input arguments are scalar. */
  for (int index=3; index<nrhs; index++)
    if(!mxIsDouble(prhs[index]) || 
      mxGetN(prhs[index])*mxGetM(prhs[index])!=1)
    {
      switch(index) {
      case 3:
        mexErrMsgTxt("Input M must be a scalar.");
      case 4:
        mexErrMsgTxt("Input xl must be a scalar.");
      case 5:
        mexErrMsgTxt("Input yl must be a scalar.");
      case 6:
        mexErrMsgTxt("Input zl must be a scalar.");
      case 7:
        mexErrMsgTxt("Input nu must be a scalar.");
      case 8:
        mexErrMsgTxt("Input mrx must be a scalar.");
      case 9:
        mexErrMsgTxt("Input mry must be a scalar.");
      case 10:
        mexErrMsgTxt("Input mrz must be a scalar.");
      case 11:
        mexErrMsgTxt("Input dt must be a scalar.");
      case 12:
        mexErrMsgTxt("Input Theta must be a scalar.");
      case 13:
        mexErrMsgTxt("Input CGTol must be a scalar.");
      case 14:
        mexErrMsgTxt("Input CGIter must be a scalar.");
      default:
        mexErrMsgTxt("Input error.");
      }
    }
  
  /* Get the matrice input arguments. */
  T0 = mxGetPr(prhs[0]);
  S1 = mxGetPr(prhs[1]);
  S2 = mxGetPr(prhs[2]);
  
  /* Get the scalar input arguments. */
  M = mxGetScalar(prhs[3]);
  xl = mxGetScalar(prhs[4]);
  yl = mxGetScalar(prhs[5]);
  zl = mxGetScalar(prhs[6]);
  nu = mxGetScalar(prhs[7]);
  mrx = mxGetScalar(prhs[8]);
  mry = mxGetScalar(prhs[9]);
  mrz = mxGetScalar(prhs[10]);
  dt = mxGetScalar(prhs[11]);
  theta = mxGetScalar(prhs[12]);
  CGTol = mxGetScalar(prhs[13]);
  CGIter = mxGetScalar(prhs[14]);
  
  /* Get the dimensions of the input matrix S1. */
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  
  // test begin
  //  for (int count=0; count<mrows*ncols; count++)
  //    mexPrintf("%g ", *(T0+count));
  // test end
  
  /* Set the output pointer to the output matrix. */
  plhs[0] = mxCreateDoubleMatrix(mrows, ncols-1, mxREAL);
  
  /* Create a C pointer to a copy of the output matrix. */
  sol2D = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. */
  ipdrops(sol2D, T0, S1, S2, M, xl, yl, zl, nu, mrx, mry, mrz, 
    dt, ncols-1, theta, CGTol, CGIter);
    
  return;
  
}

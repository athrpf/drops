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
    static int _FacePts, _SqrtFacePts;
	
  public:
    MatConnect(double DeltaT, double xl,
      double yl, double zl, int FacePts,
      double* S1, double* S2, double* S3,
      double* S4, double* S5, double* S6)
    {
      _DeltaT= DeltaT;
      _XLen= xl; _YLen= yl; _ZLen= zl;
      _SpIncrX= xl/(sqrt(FacePts)-1.0);
      _SpIncrY= yl/(sqrt(FacePts)-1.0);
      _SpIncrZ= zl/(sqrt(FacePts)-1.0);
      _FacePts= FacePts;
      _BndData[0]= S1; _BndData[1]= S2; _BndData[2]= S3;
      _BndData[3]= S4; _BndData[4]= S5; _BndData[5]= S6;
      
      _SqrtFacePts= 0;
      for (double coord= .0; coord<= xl+_SpIncrX/2.0; coord+= _SpIncrX)
        _SqrtFacePts++;
    }
    
    template<int num> static int getLexNum(const DROPS::Point2DCL& p)
    {
      int count= 0;
      if (num==0)  // Punkt liegt auf der Seite S1 oder S2
      {
        for (double zcoord= .0; zcoord<= _ZLen*p[1]-_SpIncrZ/2.0; zcoord+= _SpIncrZ)
          count+= _SqrtFacePts;
        for (double ycoord= .0; ycoord<= _YLen*p[0]+_SpIncrY/2.0; ycoord+= _SpIncrY)
          count++;
      }
      else if (num==1)  // Punkt liegt auf einer der Seiten S3 oder S4
      {
        for (double zcoord= .0; zcoord<= _ZLen*p[1]-_SpIncrZ/2.0; zcoord+= _SpIncrZ)
          count+= _SqrtFacePts;
        for (double xcoord= .0; xcoord<= _XLen*p[0]+_SpIncrX/2.0; xcoord+= _SpIncrX)
          count++;
      }
      else  // Punkt liegt auf einer der Seiten S5 oder S6
      {
        for (double ycoord= .0; ycoord<= _YLen*p[1]-_SpIncrY/2.0; ycoord+= _SpIncrY)
          count+= _SqrtFacePts;
        for (double xcoord= .0; xcoord<= _XLen*p[0]+_SpIncrX/2.0; xcoord+= _SpIncrX)
          count++;
      }
    
      return count;
    }
  
    template<int num> static double getBndVal(const DROPS::Point2DCL& p, double t)
    {
      int count= -1;  // das Feld beginnt mit Index 0
    
      // die Matrizen S1..S6 werden spaltenweise uebergeben
      for (double time= .0; time< t-_DeltaT/2.0; time+= _DeltaT)
        count+= _FacePts;
      if (num==0 || num==1)  // Seite S1 oder S2
        count+= getLexNum<0>(p);
      else if (num==2 || num==3) // Seite S3 oder S4
        count+= getLexNum<1>(p);
      else  // Seite S5 oder S6
        count+= getLexNum<2>(p);
    
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
int MatConnect::_FacePts= 0;
int MatConnect::_SqrtFacePts= 0;


namespace DROPS // for Strategy
{

template<class MGB, class Coeff>
void Strategy(InstatPoissonP1CL<MGB, Coeff>& Poisson, double nu,
  double dt, int time_steps, double cut_pos, int face_pts,
  double* T0, double* sol2D)
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
  
  // Daten fuer das PCG-Verfahren
  double tol= 1.0e-7;
  int max_iter= 500;
  
  idx.Set( 0, 1, 0, 0, 0);
  
  MultiGridCL& MG= Poisson.GetMG();
  
  // erzeuge Nummerierung zu diesem Index
  Poisson.CreateNumbering(MG.GetLastLevel(), &idx);
  
  // Vektoren mit Index idx
  b.SetIdx( &idx);
  x.SetIdx( &idx);
  cplA.SetIdx( &idx);
  cplM.SetIdx( &idx);
  
  mexPrintf("Anzahl der Unbekannten: %d\n", x.Data.size());
  
  // Steifigkeitsmatrix mit Index idx (Zeilen und Spalten)
  A.SetIdx( &idx, &idx);
  // Massematrix mit Index idx (Zeilen und Spalten)
  M.SetIdx( &idx, &idx);
  
  // stationaerer Anteil
  Poisson.SetupInstatSystem(A, M);
  
  // instationaere rechte Seite
  Poisson.SetupInstatRhs( cplA, cplM, t, b, t);
  
  // PCG-Verfahren mit SSOR-Vorkonditionierer
  SSORPcCL pc(1.0);
  PCG_SsorCL pcg_solver(pc, max_iter, tol);
  
  // Zeitdiskretisierung mit one-step-theta-scheme
  // theta=1 -> impl. Euler; theta=0.5 -> Crank-Nicholson
  InstatPoissonThetaSchemeCL<InstatPoissonP1CL<MGB, Coeff>, PCG_SsorCL>
    ThetaScheme(Poisson, pcg_solver, 0.5);
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
        for (int i=0; i<cut_pos-0.1; i++)
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
void ipdrops(double nu, double dt, double xl, double yl,
  double zl, int time_steps, int face_pts, double* T0,
  double* S1, double* S2, double* sol2D, double M)
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
    
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
    
    mexPrintf("\nRueckgabe der Daten fuer die Flaeche x=%g\n", 
      M*xl/(sqrt(face_pts)-1.0));
    mexPrintf("\nDelta t = %g", dt);
    mexPrintf("\nAnzahl der Zeitschritte = %d\n", time_steps);
    //mexPrintf("\nAnzahl der Pkte einer Seite = %d\n", face_pts);
    
    int dim= face_pts*(time_steps+1);
    double test[dim];
    for (int count=0; count<dim; count++)
      test[count]= 0.0;
    double* testp= &test[0];
    //for (int count=0; count<dim; count++)
    //  mexPrintf("%g", *(testp+count));
    //mexPrintf("\n");
    //for (int count=0; count<dim; count++)
    //  mexPrintf("%g", *(S1+count));
    //mexPrintf("\n");
    
    MatConnect MatCon(dt, xl, yl, zl, face_pts, S1, S2,
      testp, testp, testp, testp);
      
    
    const bool isneumann[6]= {true, true, true, true, true, true};
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[6]=
      { &MatConnect::getBndVal<0>, &MatConnect::getBndVal<1>,
        &MatConnect::getBndVal<2>, &MatConnect::getBndVal<3>,
        &MatConnect::getBndVal<4>, &MatConnect::getBndVal<5>};
    
    DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();
    
    int brick_div= 1;
    int pot= 2;
    int value= 9;
    
    while (value < face_pts)
    {
      brick_div++;
      pot*= 2;
      value= (pot+1)*(pot+1);
      
      MarkAll(mg);
      mg.Refine();
    }
    
    mexPrintf("Anzahl der Verfeinerungen = %d\n", brick_div);
    
    // test begin
    /*
      DROPS::Point2DCL point;
      double sp_incr= 1.0/(sqrt(face_pts)-1);
      double time= 0.0;
      for (int count=0; count<=time_steps; count++)
      {
        for (double ycount=0.0; ycount<=1.0; ycount+=sp_incr)
          for (double xcount=0.0; xcount<=1.0; xcount+=sp_incr)
          {
            point[0]= xcount;
            point[1]= ycount;
            mexPrintf("Wert an der Stelle (%g, %g) zum Zeitpkt %g: %g\n",
              point[0], point[1], time, MatCon.getBndVal<0>(point, time));
          }
        time+= dt;
      }
      
      DROPS::Point2DCL point;
      double sp_incr= 1.0/(sqrt(face_pts)-1);
      double sp_incrX= xl/(sqrt(face_pts)-1);
      double time= 0.0;
      for (int count=0; count<=time_steps; count++)
      {
        for (double ycount=0.0; ycount<=1.0; ycount+=sp_incr)
          for (double xcount=0.0; xcount<=xl; xcount+=sp_incrX)
          {
            point[0]= xcount;
            point[1]= ycount;
            mexPrintf("Wert an der Stelle (%g, %g) zum Zeitpkt %g: %g\n",
              point[0], point[1], time, MatCon.getBndVal<5>(point, time));
          }
        time+= dt;
      }
    */
    // test end
    
    //mg.SizeInfo();
    DROPS::Strategy(prob, nu, dt, time_steps, M, face_pts, T0, sol2D);
    
    return;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double nu, dt, xl, yl, zl, M;
  double *T0, *S1, *S2, *sol2D;
  int mrows, ncols;
  
  /* Check for proper number of arguments. */
  if(nrhs!=9)
    mexErrMsgTxt("(nu, dt, xl, yl, zl, T0, S1, S2, M) as input required.");
  if(nlhs!=1)
    mexErrMsgTxt("Solution on 2D-area as output required.");
    
  /* Check to make sure the first input argument is a scalar. */
  if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
    mxIsEmpty(prhs[0]) || mxIsComplex(prhs[0]) ||
    mxGetN(prhs[0])*mxGetM(prhs[0])!=1 )
  {
    mexErrMsgTxt("Input nu must be a scalar.");
  }
  
  /* Check to make sure the second input argument is a scalar. */
  if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
    mxIsEmpty(prhs[1]) || mxIsComplex(prhs[1]) ||
    mxGetN(prhs[1])*mxGetM(prhs[1])!=1 )
  {
    mexErrMsgTxt("Input dt must be a scalar.");
  }
  
  /* Check to make sure the third input argument is a scalar. */
  if( !mxIsNumeric(prhs[2]) || !mxIsDouble(prhs[2]) ||
    mxIsEmpty(prhs[2]) || mxIsComplex(prhs[2]) ||
    mxGetN(prhs[2])*mxGetM(prhs[2])!=1 )
  {
    mexErrMsgTxt("Input xl must be a scalar.");
  }
  
  /* Check to make sure the fourth input argument is a scalar. */
  if( !mxIsNumeric(prhs[3]) || !mxIsDouble(prhs[3]) ||
    mxIsEmpty(prhs[3]) || mxIsComplex(prhs[3]) ||
    mxGetN(prhs[3])*mxGetM(prhs[3])!=1 )
  {
    mexErrMsgTxt("Input yl must be a scalar.");
  }
  
  /* Check to make sure the fifth input argument is a scalar. */
  if( !mxIsNumeric(prhs[4]) || !mxIsDouble(prhs[4]) ||
    mxIsEmpty(prhs[4]) || mxIsComplex(prhs[4]) ||
    mxGetN(prhs[4])*mxGetM(prhs[4])!=1 )
  {
    mexErrMsgTxt("Input zl must be a scalar.");
  }
  
  /* Check to make sure the last input argument is a scalar. */
  if( !mxIsNumeric(prhs[8]) || !mxIsDouble(prhs[8]) ||
    mxIsEmpty(prhs[8]) || mxIsComplex(prhs[8]) ||
    mxGetN(prhs[8])*mxGetM(prhs[8])!=1 )
  {
    mexErrMsgTxt("Input M must be a scalar.");
  }
  
  /* Get the scalar input nu. */
  nu = mxGetScalar(prhs[0]);
  
  /* Get the scalar input dt. */
  dt = mxGetScalar(prhs[1]);
  
  /* Get the scalar input xl. */
  xl = mxGetScalar(prhs[2]);
  
  /* Get the scalar input yl. */
  yl = mxGetScalar(prhs[3]);
  
  /* Get the scalar input zl. */
  zl = mxGetScalar(prhs[4]);
  
  /* Create a pointer to the input matrices T0,S1,S2. */
  T0 = mxGetPr(prhs[5]);
  S1 = mxGetPr(prhs[6]);
  S2 = mxGetPr(prhs[7]);
  
  /* Get the scalar input M. */
  M = mxGetScalar(prhs[8]);
  
  /* Get the dimensions of the input matrix S1. */
  mrows = mxGetM(prhs[6]);
  ncols = mxGetN(prhs[6]);
  
  // test begin
  
  //  for (int count=0; count<mrows*ncols; count++)
  //    mexPrintf("%g ", *(T0+count));
  
  // test end
  
  /* Set the output pointer to the output matrix. */
  plhs[0] = mxCreateDoubleMatrix(mrows, ncols-1, mxREAL);
  
  /* Create a C pointer to a copy of the output matrix. */
  sol2D = mxGetPr(plhs[0]);
  
  /* Call the C subroutine. */
  ipdrops(nu, dt, xl, yl, zl, ncols-1, mrows, T0, S1, S2, sol2D, M);
  
  return;
  
}

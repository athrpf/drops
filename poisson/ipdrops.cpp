//**************************************************************************
// File:    ipdrops.cpp                                                    *
// Content: test program for the instat. poisson-problem                   *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, Marcus Soemers      *
//          IGPM RWTH Aachen                                               *
// Version: 0.1                                                            *
// History: begin - Nov, 19 2002                                           *
//**************************************************************************

#include "poisson/instatpoisson.h"
#include "num/poissonsolver.h"
#include "poisson/integrTime.h"


// du/dt - laplace u + q*u = f

class PoissonCoeffCL
{
  public:
    static double q(const DROPS::Point3DCL&) { return 0.0; }
    static double f(const DROPS::Point3DCL& , double )
      { return 0.0; }
    //static double f(const DROPS::Point3DCL& p, double t)
    //  { return (-2.0*exp(t)*exp(p[0]+p[1]+p[2])); }
};


inline double Lsg(const DROPS::Point3DCL& , double )
  { return (0.0); }
//inline double Lsg(const DROPS::Point3DCL& p, double t)
//  { return (exp(t)*exp(p[0]+p[1]+p[2])); }


// boundary functions (neumann, dirichlet type)
// used for BndSegCL-object of a UnitCube


/*  
inline double neu_val1(const DROPS::Point2DCL& p, double t)
  { return (-exp(t)*exp(p[0]+p[1])); }
inline double neu_val2(const DROPS::Point2DCL& p, double t)
  { return (exp(t)*exp(1.0+p[0]+p[1])); }
*/
inline double neu_val1(const DROPS::Point2DCL& p, double t)
  { return (-4.0); }
inline double neu_val2(const DROPS::Point2DCL& p, double t)
  { return (4.0); }
inline double neu_val3(const DROPS::Point2DCL& p, double t)
  { return (0.0); }
  
  
inline double dir_val1(const DROPS::Point2DCL& p, double t)
  { return (exp(t)*exp(p[0]+p[1])); }
inline double dir_val2(const DROPS::Point2DCL& p, double t)
  { return (exp(t)*exp(1.0+p[0]+p[1])); }
  
  

namespace DROPS // for Strategy
{

template<class MGB, class Coeff>
void Strategy(InstatPoissonP1CL<MGB, Coeff>& Poisson, double dt,
  double time_steps)
{	
  typedef InstatPoissonP1CL<MGB,Coeff> MyPoissonCL;
  
  IdxDescCL& idx= Poisson.idx;
  VecDescCL& x= Poisson.x;
  VecDescCL& b= Poisson.b;
  MatDescCL& A= Poisson.A;
  MatDescCL& M= Poisson.M;
  
  VecDescCL cplA;
  VecDescCL cplM;
	
  // Zeitschrittweite, aktueller Zeitpunkt
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
  
  std::cerr << "Anzahl der Unbekannten: " <<  x.Data.size() << std::endl;
    
  // Steifigkeitsmatrix mit Index idx (Zeilen und Spalten)
  A.SetIdx( &idx, &idx);        
  // Massematrix mit Index idx (Zeilen und Spalten)
  M.SetIdx( &idx, &idx);
  
  // stationaerer Anteil
  Poisson.SetupInstatSystem(A, M);
  
  // instationaere rechte Seite 
  Poisson.SetupInstatRhs( cplA, cplM, t, b, t);
   
  SsorPcCL<VectorCL, double> pc(1.0);
  PCG_T pcg_solver(tol, max_iter, pc);
  InstatPoissonThetaSchemeCL<InstatPoissonP1CL<MGB, Coeff> ,PCG_T> 
    ThetaScheme(Poisson, pcg_solver, 0.5);
  
  ThetaScheme.SetTimeStep(dt);
  
  scalar_instat_fun_ptr exact_sol = &Lsg;
  
  
  // ****** Startwert
  
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
    
  for (MultiGridCL::TriangVertexIteratorCL sit=MG.GetTriangVertexBegin(lvl), 
    send=MG.GetTriangVertexEnd(lvl); sit != send; ++sit)
  { 
    if (sit->Unknowns.Exist())
    {
      IdxT i= sit->Unknowns(indx)[0];
      pt= sit->GetCoord();
    
      help= std::make_pair(pt[2], pt[1]);
      key= std::make_pair(pt[0], help);
      x.Data[i]= Lsg(pt, 0.0);
      nmap[key]= &(x.Data[i]);
    }
  }
  
  
  // Ausgabe Startwert
  /*
  for (ci p= nmap.begin(); p!= nmap.end(); p++)
  {
    std::cerr << *(p->second) << "\n";
  }
  */
  
  
  // ****** Ende Startwert
  
    
  for (int step=1;step<=time_steps;step++)
  {
    ThetaScheme.DoStep(x);
    std::cerr << "t= " << Poisson.t << std::endl;
    std::cerr << "Iterationen: " << pcg_solver.GetIter() 
      << "    Norm des Residuums: " << pcg_solver.GetResid() << std::endl;
    Poisson.CheckSolution(x, exact_sol, Poisson.t);
  }
  
  A.Reset();
  b.Reset();
  
  
  // Ausgabe Loesung   
  /*
  for (ci p= nmap.begin(); p!= nmap.end(); p++)
  {
    std::cerr << *(p->second) << "\n";
  }
  */
  
}

} // end of namespace DROPS



int main()
{
  try
  {
    DROPS::Point3DCL null(0.0);
    DROPS::Point3DCL e1(0.0), e2(0.0), e3(0.0);
    e1[0]= 1.0;
    e2[1]= e3[2]= 1.0;

    typedef DROPS::InstatPoissonP1CL<DROPS::BrickBuilderCL, PoissonCoeffCL> 
      InstatPoissonOnBrickCL;
    typedef InstatPoissonOnBrickCL MyPoissonCL;
    
    DROPS::BrickBuilderCL brick(null, e1, e2, e3, 2, 2, 2);
  
    double dt;
    int time_steps, brick_div;
  
    std::cerr << "\nDelta t = "; std::cin >> dt;
    std::cerr << "\nAnzahl der Zeitschritte = "; std::cin >> time_steps;
    std::cerr << "\nAnzahl der Verfeinerungen = "; std::cin >> brick_div;

    // Neumann boundary conditions   
    const bool isneumann[6]= 
      { true, true, true, true, true, true };
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[24]=
      { &neu_val1, &neu_val2, &neu_val3, &neu_val3, &neu_val3, &neu_val3 };
      /*
    const bool isneumann[6]= 
      { true, true, true, true, true, true };
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[24]=
      { &neu_val1, &neu_val2, &neu_val1, &neu_val2, &neu_val1, &neu_val2 };
  */
  
    /*
    // Dirichlet boundary conditions 
    const bool isneumann[6]= 
      { false, false, false, false, false, false };
    const DROPS::InstatPoissonBndDataCL::bnd_val_fun bnd_fun[24]=
      { &dir_val1, &dir_val2, &dir_val1, &dir_val2, &dir_val1, &dir_val2 };
    */
      
    DROPS::InstatPoissonBndDataCL bdata(6, isneumann, bnd_fun);
    MyPoissonCL prob(brick, PoissonCoeffCL(), bdata);
    DROPS::MultiGridCL& mg = prob.GetMG();
    
    for (int count=1; count<brick_div; count++)
    {
      MarkAll(mg);
      mg.Refine();
    } 
    
    mg.SizeInfo(std::cerr);
    DROPS::Strategy(prob, dt, time_steps);
   
   
    return 0;
  }
  catch (DROPS::DROPSErrCL err) { err.handle(); }
  
}




//**************************************************************************
// File:    params.h                                                       *
// Content: parameters for two-phase flow problems                         *
// Author:  Sven Gross, Joerg Peters, Volker Reichelt, IGPM RWTH Aachen    *
//**************************************************************************

#ifndef DROPS_LSET_PARAMS_H
#define DROPS_LSET_PARAMS_H

#include "misc/params.h"

namespace DROPS
{

class ParamMesszelleCL: public ParamBaseCL
{
  private:
    void RegisterParams();
    
  public:
    int    StokesMethod;                        // Loeserklasse fuer die Stokesprobleme
    double inner_tol, outer_tol, 		// Parameter der Loeser
           lset_tol, lset_SD;			// fuer Flow & Levelset
    int    inner_iter, outer_iter, lset_iter;	
    int    FPsteps;				// Kopplung Levelset/Flow: Anzahl Fixpunkt-Iterationen

    double dt;					// Zeitschrittweite
    int    num_steps;				// Anzahl Zeitschritte
    double theta, lset_theta;			// 0=FwdEuler, 1=BwdEuler, 0.5=CN

    double sigma,				// Oberflaechenspannung
           CurvDiff,				// num. Glaettung Kruemmungstermberechnung
           rhoD, rhoF, muD, muF, 	 	// Stoffdaten: Dichte/Viskositaet
           sm_eps; 				// Glaettungszone fuer Dichte-/Viskositaetssprung

    Point3DCL g;				// Schwerkraft
    double    Radius;				// Radius und 
    Point3DCL Mitte;				// Position des Tropfens
    int       num_dropref,			// zusaetzliche Tropfenverfeinerung
              VolCorr,				// Volumenkorrektur (0=false)
              IniCond;				// Anfangsbedingung (0=Null, 1/2= stat. flow mit/ohne Tropfen)
              
    int    ref_flevel, ref_freq;		// Parameter fuer
    double ref_width;				// adaptive Verfeinerung
    
    double Anstroem, 				// max. Einstroemgeschwindigkeit (Parabelprofil)
           r_inlet;				// Radius am Einlass der Messzelle
    int    flow_dir;				// Stroemungsrichtung (x/y/z = 0/1/2)

    int    RepFreq, RepSteps, RepMethod;	// Parameter fuer
    double RepTau, RepDiff;  			// Reparametrisierung

    string EnsCase,				// Ensight Case, 
           EnsDir,				// lok.Verzeichnis, in das die geom/vec/scl-files abgelegt werden
           IniData,
           meshfile;				// Meshfile (von GAMBIT im FLUENT/UNS-Format)

    ParamMesszelleCL()                        { RegisterParams(); }
    ParamMesszelleCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

class ParamMesszelleNsCL: public ParamMesszelleCL
{
  private:
    void RegisterParams();
    
  public:
    int    scheme;				// 0=Baensch, 1=theta-scheme
    double nonlinear,				// Anteil Nichtlinearitaet
           stat_nonlinear, stat_theta;		// stat. Rechnung fuer Anfangswerte
  
    ParamMesszelleNsCL()
      : ParamMesszelleCL() { RegisterParams(); }
    ParamMesszelleNsCL( const string& filename) 
      : ParamMesszelleCL() { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

class ParamFilmCL: public ParamBaseCL
{ // x = Filmnormal, z = Ablaufrichtung
  private:
    void RegisterParams();
    
  public:
    double inner_tol, outer_tol, 		// Parameter der Loeser
           lset_tol, lset_SD;			// fuer Flow & Levelset
    int    inner_iter, outer_iter, lset_iter;	
    int    FPsteps;				// Kopplung Levelset/Flow: Anzahl Fixpunkt-Iterationen

    double dt;					// Zeitschrittweite
    int    num_steps;				// Anzahl Zeitschritte
    double theta, lset_theta;			// 0=FwdEuler, 1=BwdEuler, 0.5=CN

    double sigma,				// Oberflaechenspannung
           CurvDiff,				// num. Glaettung Kruemmungstermberechnung
           rhoF, rhoG, muF, muG, 	 	// Stoffdaten: Dichte/Viskositaet
           sm_eps, 				// Glaettungszone fuer Dichte-/Viskositaetssprung
           PumpAmpl, PumpFreq;			// Frequenz und Amplitude der Anregung

    Point3DCL g;				// Schwerkraft
    double    Filmdicke;			// Filmdicke
    Point3DCL mesh_res,				// Gitteraufloesung und
              mesh_size;			// Gittergroesse in x-/y-/z-Richtung
    int       num_ref,				// zusaetzliche Verfeinerung
              VolCorr,				// Volumenkorrektur (0=false)
              IniCond;				// Anfangsbedingung (0=Null, 1= stat. flow, -1= read from file )
              
    int    ref_flevel, ref_freq;		// Parameter fuer
    double ref_width;				// adaptive Verfeinerung
    
    int    RepFreq, RepMethod;			// Parameter fuer Reparametrisierung

    string EnsCase,				// Ensight Case, 
           EnsDir,				// lok.Verzeichnis, in das die geom/vec/scl-files abgelegt werden
           IniData,
           BndCond;

    ParamFilmCL()                        { RegisterParams(); }
    ParamFilmCL( const string& filename) { RegisterParams(); std::ifstream file(filename.c_str()); rp_.ReadParams( file); }
};

} // end of namespace DROPS

#endif
    
    
    
    

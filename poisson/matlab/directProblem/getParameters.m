% Diese Funktion liefert saemtliche Parameter, die fuer den Aufruf von DROPS notwendig sind
%-------------------------------------------------------------------------------------------

function [C_data,func,L_data]=getParameters()

% Geometrieparameter
xl= 2.0;
yl= 1.0;
zl= 1.0;
gl= [xl,yl,zl];

% Skalierung
% nu=lambda/(rho*c)
nu= 1.0;

% Quellterm
% 0->ohne, 1->mit
Flag= 0;

% Zeitdiskretisierung
% 0 expl. Euler, 0.5 Crank-Nicholson, 1 impl. Euler
theta= 1.0;

% Zeitschrittweite
dt= 0.1;

% Anzahl Zeitschritte
ndt= 10;

% exakter Waermestrom
func= 8;

% ni* =   Anzahl der Intervalle pro Koordinatenrichtung
%         nach uniformer Verfeinerung
nint= 4;
nix= nint;
niy= nint;
niz= nint;

ni= [nix,niy,niz];
npyz= (niy+1)*(niz+1);

% Anfangs- und Randdaten
T0= 2;
qc= 0.1;

% Loeser fuer die DPs
% 0 PCG, 1 MG, 2 Block-GS, 3 SSOR
Solver= 0;

% allgemeine Daten fuer DROPS zusammenfassen
C_data= [gl,ni,dt,ndt,nu,Flag,theta,T0,qc,Solver];

% allgemeines Abbruchkriterium
SolverMaxIter= 1000;

%-------------------- local data ---------------------

% Schnittflaeche durch das Gebiet
L_plane= 1;

L_BndRef= 0;

% L_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung
%           vor der Verfeinerung
L_nixs= nint;
L_niys= nint;
L_nizs= nint;

L_nis= [L_nixs,L_niys,L_nizs];

L_tol= 1e-5;
L_iter= SolverMaxIter;

L_data= [L_nis,L_BndRef,L_tol,L_iter,L_plane];

return

% Diese Funktion liefert saemtliche Parameter, die fuer den Aufruf von DROPS notwendig sind
%-------------------------------------------------------------------------------------------

function [C_data,func,niter_max,epsT,M_data,D_data,A_data,S_data]=getParameters()


% Geometrieparameter
xl= 0.025;
yl= 19.5;
zl= 39;
gl= [xl,yl,zl];

% Skalierung
% nu=lambda/(rho*c)
nu= 6.303096;

% Quellterm
% 0->ohne, 1->mit
Flag= 0;

% Zeitdiskretisierung
% 0 expl. Euler, 0.5 Crank-Nicholson, 1 impl. Euler
theta= 0.5;

% Zeitschrittweite
dt= 2e-3;

% Anzahl Zeitschritte
ndt= 100;

% exakter Waermestrom
func= 4;

% ni* =   Anzahl der Intervalle pro Koordinatenrichtung
%         nach uniformer Verfeinerung
nix= 2;
niy= 4;
niz= 20;
% nix= 4;
% niy= 40;
% niz= 80;
% nix= 4;
% niy= 96;
% niz= 196;
ni= [nix,niy,niz];
npyz= (niy+1)*(niz+1);

% Anfangs- und Randdaten
T0= 36.8097; %mittlere Fluidtemperatur
qh= -278.652e-3;

% Loeser fuer die DPs
% 0 PCG, 1 MG, 2 Block-GS, 3 SSOR
Solver= 0;

% allgemeine Daten fuer DROPS zusammenfassen
C_data= [gl,ni,dt,ndt,nu,Flag,theta,T0,qh,Solver];

% maximale Anzahl der Optimierungsschritte
niter_max= 20;

% Abbruch ueber Residuum
epsT= 1e-14;

% allgemeines Abbruchkriterium
SolverMaxIter= 200;

%-------------------- measurement data ---------------------

% Schnittflaeche durch das Gebiet
M_plane= 0;

M_BndRef= 0;

% M_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung
%           vor der Verfeinerung
M_nixs= 2;
M_niys= 4;
M_nizs= 20;
% M_nixs= 4;
% M_niys= 40;
% M_nizs= 80;
% M_nixs= 1;
% M_niys= 24;
% M_nizs= 49;
M_nis= [M_nixs,M_niys,M_nizs];

M_tol= 1e-6;
M_iter= SolverMaxIter;

M_data= [M_nis,M_BndRef,M_tol,M_iter,M_plane];

%-------------------- direct problem -----------------------

% Schnittflaeche durch das Gebiet
D_plane= 0;

D_BndRef= 0;

% D_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung
%           vor der Verfeinerung
D_nixs= 2;
D_niys= 4;
D_nizs= 20;
% D_nixs= 4;
% D_niys= 40;
% D_nizs= 80;
% D_nixs= 1;
% D_niys= 24;
% D_nizs= 49;
D_nis= [D_nixs,D_niys,D_nizs];

D_tol= 1e-1;
D_iter= SolverMaxIter;

D_data= [D_nis,D_BndRef,D_tol,D_iter,D_plane];

%-------------------- adjoint problem ----------------------

% Schnittflaeche durch das Gebiet
A_plane= xl;

A_BndRef= 0;

% A_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung
%           vor der Verfeinerung
A_nixs= 2;
A_niys= 4;
A_nizs= 20;
% A_nixs= 4;
% A_niys= 40;
% A_nizs= 80;
% A_nixs= 1;
% A_niys= 24;
% A_nizs= 49;
A_nis= [A_nixs,A_niys,A_nizs];

A_tol= 1e-1;
A_iter= SolverMaxIter;

A_data= [A_nis,A_BndRef,A_tol,A_iter,A_plane];

%------------------ sensitivity problem --------------------

% Schnittflaeche durch das Gebiet
S_plane= 0;

S_BndRef= 0;

% S_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung
%           vor der Verfeinerung
S_nixs= 2;
S_niys= 4;
S_nizs= 20;
% S_nixs= 4;
% S_niys= 40;
% S_nizs= 80;
% S_nixs= 1;
% S_niys= 24;
% S_nizs= 49;
S_nis= [S_nixs,S_niys,S_nizs];

S_tol= 1e-1;
S_iter= SolverMaxIter;

S_data= [S_nis,S_BndRef,S_tol,S_iter,S_plane];

return

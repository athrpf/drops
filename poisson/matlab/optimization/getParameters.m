% Diese Funktion liefert saemtliche Parameter, die fuer den Aufruf von DROPS notwendig sind 
%-------------------------------------------------------------------------------------------

function [C_data,func,niter_max,epsT,M_data,D_data,A_data,S_data]=getParameters()

% Geometrieparameter 
xl= 10;
yl= 40;
zl= 100;
gl= [xl,yl,zl];

% Skalierung 
nu= 100;

% Quellterm
Flag= 0;

% Zeitdiskretisierung 
theta= 1.0;

% Zeitschrittweite
dt= 0.01;

% Anzahl Zeitschritte
ndt= 200;

% exakter Waermestrom 
func= 3;

% ni* =   Anzahl der Intervalle pro Koordinatenrichtung 
%         nach uniformer Verfeinerung
nix= 4;
niy= 4;
niz= 4;
ni= [nix,niy,niz];

% allgemeine Daten fuer DROPS zusammenfassen
C_data= [gl,ni,dt,ndt,nu,Flag,theta];

% maximale Anzahl der Optimierungsschritte
niter_max= 5;

% Abbruch ueber Residuum
epsT= 1e-16;


%-------------------- measurement data ---------------------

% Schnittflaeche durch das Gebiet
M_plane= 0;

M_BndRef= 0;

% M_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung 
M_nixs= 4;
M_niys= 4;
M_nizs= 4;
M_nis= [M_nixs,M_niys,M_nizs];

M_cgtol= 10^(-12);
M_cgiter= 500;

M_data= [M_nis,M_BndRef,M_cgtol,M_cgiter,M_plane];

%-------------------- direct problem -----------------------

% Schnittflaeche durch das Gebiet
D_plane= 0;

D_BndRef= 0;

% D_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung
D_nixs= 4;
D_niys= 4;
D_nizs= 4;
D_nis= [D_nixs,D_niys,D_nizs];

D_cgtol= 10^(-6);
D_cgiter= 500;

D_data= [D_nis,D_BndRef,D_cgtol,D_cgiter,D_plane];

%-------------------- adjoint problem ----------------------

% Schnittflaeche durch das Gebiet
A_plane= xl;

A_BndRef= 0;

% A_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung
A_nixs= 4;
A_niys= 4;
A_nizs= 4;
A_nis= [A_nixs,A_niys,A_nizs];

A_cgtol= 10^(-6);
A_cgiter= 500;

A_data= [A_nis,A_BndRef,A_cgtol,A_cgiter,A_plane];

%------------------ sensitivity problem --------------------

% Schnittflaeche durch das Gebiet
S_plane= 0;

S_BndRef= 0;

% S_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung
S_nixs= 4;
S_niys= 4;
S_nizs= 4;
S_nis= [S_nixs,S_niys,S_nizs];

S_cgtol= 10^(-6);
S_cgiter= 500;

S_data= [S_nis,S_BndRef,S_cgtol,S_cgiter,S_plane];

return

% Diese Funktion liefert saemtliche Parameter, die fuer den Aufruf von DROPS notwendig sind 
%-------------------------------------------------------------------------------------------

function [C_data,func,niter_max,epsT,M_data,D_data,A_data,S_data]=getParameters()

% Geometrieparameter 
xl= 0.025;
yl= 15;
zl= 30;
gl= [xl,yl,zl];

% Skalierung 
nu= 6.303096;

% Quellterm
Flag= 0;

% Zeitdiskretisierung 
theta= 1.0;

% Zeitschrittweite
dt= 0.000666666;

% Anzahl Zeitschritte
ndt= 190;

% exakter Waermestrom 
func= 3;

% ni* =   Anzahl der Intervalle pro Koordinatenrichtung 
%         nach uniformer Verfeinerung
nix= 4;
niy= 63;
niz= 127;
ni= [nix,niy,niz];

% allgemeine Daten fuer DROPS zusammenfassen
C_data= [gl,ni,dt,ndt,nu,Flag,theta];

% maximale Anzahl der Optimierungsschritte
niter_max= 10;

% Abbruch ueber Residuum
epsT= 1e-16;


%-------------------- measurement data ---------------------

% Schnittflaeche durch das Gebiet
M_plane= 0;

M_BndRef= 0;

% M_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung 
M_nixs= 4;
M_niys= 63;
M_nizs= 127;
M_nis= [M_nixs,M_niys,M_nizs];

M_cgtol= 1e-16;
M_cgiter= 500;

M_data= [M_nis,M_BndRef,M_cgtol,M_cgiter,M_plane];

%-------------------- direct problem -----------------------

% Schnittflaeche durch das Gebiet
D_plane= 0;

D_BndRef= 0;

% D_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung
D_nixs= 4;
D_niys= 63;
D_nizs= 127;
D_nis= [D_nixs,D_niys,D_nizs];

D_cgtol= 1e-14;
D_cgiter= 500;

D_data= [D_nis,D_BndRef,D_cgtol,D_cgiter,D_plane];

%-------------------- adjoint problem ----------------------

% Schnittflaeche durch das Gebiet
A_plane= xl;

A_BndRef= 0;

% A_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung
A_nixs= 4;
A_niys= 63;
A_nizs= 127;
A_nis= [A_nixs,A_niys,A_nizs];

A_cgtol= 1e-14;
A_cgiter= 500;

A_data= [A_nis,A_BndRef,A_cgtol,A_cgiter,A_plane];

%------------------ sensitivity problem --------------------

% Schnittflaeche durch das Gebiet
S_plane= 0;

S_BndRef= 0;

% S_ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%           vor der Verfeinerung
S_nixs= 4;
S_niys= 63;
S_nizs= 127;
S_nis= [S_nixs,S_niys,S_nizs];

S_cgtol= 1e-14;
S_cgiter= 500;

S_data= [S_nis,S_BndRef,S_cgtol,S_cgiter,S_plane];

return

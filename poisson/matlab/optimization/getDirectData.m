function [nu,xl,yl,zl,M,theta,cgtol,cgiter,Flag,BndRef,nis,ni,func]=getDirectData()

% Parameter
nu= 100;
xl= 10;
yl= 40;
zl= 100;

BndRef= 2;

% Parameter waehlt den exakten Waermestrom aus
% func=3: sinusfoermig, stationaer 
func= 3;

% ni*s =  Anzahl der Intervalle pro Koordinatenrichtung 
%         vor der Verfeinerung 
nixs= 2;
niys= 2;
nizs= 2;
nis= [nixs,niys,nizs];

% ni* =   Anzahl der Intervalle pro Koordinatenrichtung 
%         nach uniformer Verfeinerung
nix= nixs*(2^BndRef);
niy= niys*(2^BndRef);
niz= nizs*(2^BndRef);
ni= [nix,niy,niz];

M= 0;
theta= 1.0;
cgtol= 10^(-3);
cgiter= 500;
Flag= 0;

return

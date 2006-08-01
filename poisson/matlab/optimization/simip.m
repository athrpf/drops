function [Iter,Ti]=simip(qc,C_data,L_data)

% ------------------------------------------------------------------------
% PROBLEM: LOESUNG DER GLEICHUNG
%           dT/dt - nu * laplace(T) = 0
% T: Temperatur
% t: Zeit
% nu : konstanter Parameter
% ------------------------------------------------------------------------
%
% GEOMETRIE und RAUMDISKRETISIERUNG
% Laenge in y- bzw. z-Richtung = 1
% xl= Laenge des Quaders in x-Richtung
% mr = Anzahl der Verfeinerungen
% npc = Anzahl der Punkte pro Koordinatenrichtung
% nfp= Anzahl der Punkte pro Seitenflaeche
% ------------------------------------------------------------------------
%
% ANFANGSTEMPERATUR
% lexikographische Nummerierung der Gitterpunkte
% von x=0 bis x=1
%
% ------------------------------------------------------------------------
% RANDBEDINGUNGEN
% Null auf allen Seiten ausser x=0 und x=1
% qh: Newman auf Wuerfelseite x=0
% qc: Newman auf Wuerfelseite x=1
% lexikographische Nummerierung der Gitterpunkte
%
% ------------------------------------------------------------------------
% ZEITDISKRETISIERUNG
% dt: Zeitschrittweite
% ndt: Anzahl der Zeitschritte
%
% ------------------------------------------------------------------------
% MESSPUNKTE
% M: Parameter, der die Flaeche spezifiziert an der die
% Messungen durchgefuehrt werden
% (Schnitt durch den Wuerfel an der Stelle x=M)
%--------------------------------------------------------------------------

% set common data
xl= C_data(1);
yl= C_data(2);
zl= C_data(3);
nix= C_data(4);
niy= C_data(5);
niz= C_data(6);
dt= C_data(7);
ndt= C_data(8);
nu= C_data(9);
Flag= C_data(10);
theta= C_data(11);
Solver= C_data(14);

npx= nix+1;
npyz= (niy+1)*(niz+1);

% set local data
nis= [L_data(1),L_data(2),L_data(3)];
BndRef= L_data(4);
cgtol= L_data(5);
cgiter= L_data(6);
M= L_data(7);

T0= zeros(npyz, npx);
qh= zeros(npyz, ndt+1);

% Loesung mit DROPS
work_dir= cd;
cd ..;
[Val1,T]= ipdrops(T0, qh, qc, M, xl, yl, zl, nu, nis(1), nis(2), nis(3), dt, theta, cgtol, cgiter, Flag, BndRef);
cd( work_dir);

Plane_Pos= round(M*nix/xl);
Ti= [T0(:,Plane_Pos+1) T];
Iter= Val1;

%figure(2);plotsol(npc,ndt,sol)

return

function Td=simdp()


% ------------------------------------------------------------------------
% PROBLEM: LOESUNG DER GLEICHUNG 
%           dT/dt - nu * laplace(T) = 0
% T: Temperatur
% t: Zeit
% nu : konstanter Parameter
% ------------------------------------------------------------------------
%
% GEOMETRIE und RAUMDISKRETISIERUNG
% xl,yl und zl = Laenge des Quaders in x-, y- und z-Richtung
% npx = Anzahl der Punkte in die x-Koordinatenrichtung
% npyz= Anzahl der Punkte auf yz-Seitenflaeche 
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

% gemeinsame Daten der direkten Probleme holen
[nu,xl,yl,zl,M,theta,cgtol,cgiter,Flag,BndRef,nis,ni,func]=getDirectData;

npx= ni(1)+1;
npyz= (ni(2)+1)*(ni(3)+1);
%
dt= 0.01;
ndt= 100;
%time=[0:dt:ndt*dt];
%
T0= 20*ones(npyz,npx);
%
qh= 2*ones(npyz, ndt+1);
qc=zeros(npyz, ndt+1);

% Loesung mit DROPS
cd ..;
T= ipdrops(T0, qh, qc, M, xl, yl, zl, nu, nis(1), nis(2), nis(3), dt, theta, cgtol, cgiter, Flag, BndRef);
cd optimization;
Td=[T0(:,M+1) T];
%figure(2);plotsol(npc,ndt,sol)

function Ti=simipsens(qc)


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

% gemeinsame Daten der direkten Probleme holen
[nu,xl,yl,zl,M,theta,cgtol,cgiter,Flag,BndRef,nis,ni,func]=getSensitivityData;

npx= (ni(1)+1);
npyz= (ni(2)+1)*(ni(3)+1);

dt= 0.01;
ndt= 100;
time=[0:dt:ndt*dt];
%
%npc= 2^mr+1;
%nfp= (npc)^2;
%
%T0= 20*ones(nfp, npc);
%qh= 2*ones(nfp, ndt+1);
T0= zeros(npyz, npx);
qh= zeros(npyz, ndt+1);
% qcpv=qcp*ones(1,npc);
% qcs=[];    
% for i=1:npc
%     qcs=[qcs;qcpv(npc-i+1,:)'];
% end
% qc=[];
% for i=1:ndt+1
%     qc(:,i)=qcs;
% end
%qc= qcp*ones(nfp, ndt+1);

% Loesung mit DROPS
cd ..;
T= ipdrops(T0, qh, qc, M, xl, yl, zl, nu, nis(1), nis(2), nis(3), dt, theta, cgtol, cgiter, Flag, BndRef);
cd optimization;
Ti=[T0(:,M+1) T];

%figure(2);plotsol(npc,ndt,sol)

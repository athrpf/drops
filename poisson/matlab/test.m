clear;

%------------------------------ Parameter ---------------------------------- 

% Diskretisierung der Gleichung 
% du/dt - nu * laplace(u) = 0
% setze Parameter nu=const
nu= 10.0;

% Zeitschrittweite
dt= 0.001;

% Laenge des Quaders in x-, y- und z-Richtung
xl= 1; 
yl= 1;
zl= 1;

% Anzahl der Verfeinerungen
mr= 2;

% Anzahl der Zeitschritte
ndt= 100;

% Parameter, der die Flaeche spezifiziert an der die 
% Messungen durchgefuehrt werden
% (Schnitt durch den Wuerfel an der Stelle x=M)
M= 0;



% ----------------------------- Eingabedaten -------------------------------

% Anzahl der Punkte pro Koordinatenrichtung
npc= 2^mr+1;

% Anzahl der Punkte pro Seitenflaeche 
nfp= (npc)^2;

% Anfangstemperatur
% lexikographische Nummerierung der Gitterpunkte
% von x=0 bis x=1
T0= zeros(nfp, npc);
%T0= [1 10 19
%     2 11 20
%     3 12 21
%     4 13 22
%     5 14 23
%     6 15 24
%     7 16 25
%     8 17 26
%     9 18 27];

% Daten auf den beiden gegenueberliegenden
% Wuerfelseiten (x=0 & x=1)
% lexikographische Nummerierung der Gitterpunkte
% fuer x=0 bzw. x=1
S1= -4*ones(nfp, ndt+1);
S2= 4*ones(nfp, ndt+1);

%S1= [1 10 19 28 37
%     2 11 20 29 38
%     3 12 21 30 39
%     4 13 22 31 40
%     5 14 23 32 41
%     6 15 24 33 42
%     7 16 25 34 43
%     8 17 26 35 44
%     9 18 27 36 45]; 
     
%S1= [1 10 19 28
%    2 11 20 29
%    3 12 21 30
%    4 13 22 31
%    5 14 23 32
%    6 15 24 33
%    7 16 25 34
%    8 17 26 35
%    9 18 27 36]; 

%T0= zeros(9, 81);
%S1= 4*ones(81, 11);
%S2= -4*ones(81, 11);


% Loesung mit DROPS
sol= ipdrops(nu, dt, xl, yl, zl, T0, S1, S2, M);


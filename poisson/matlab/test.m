clear;

%------------------------------ Parameter ---------------------------------- 

% Diskretisierung der Gleichung 
% du/dt - nu * laplace(u) = 0
% setze Parameter nu=const
nu= 10.0;

% Zeitschrittweite
dt= 0.001;

% Laenge des Quaders in x-, y- und z-Richtung
xl= 1.0; 
yl= 1.0;
zl= 1.0;

% Anzahl der Teilintervalle pro Koordinatenrichtung 
mrx= 4;
mry= 4;
mrz= 4;

% Anzahl der Zeitschritte
ndt= 100;

% Zeitdiskretisierung ueber one-step-theta-scheme
% theta= 1.0 -> impliziter Euler
% theta= 0.5 -> Crank-Nicholson 
% theta= 0.0 -> expliziter Euler
theta= 0.5;

% Toleranz fuer das CG-Verfahren
cgtol= 10^(-4);
% maximale Anzahl der CG-Iterationen
cgiter= 500;



% Parameter, der die Flaeche spezifiziert an der die 
% Messungen durchgefuehrt werden
% (Schnitt durch den Wuerfel an der Stelle x=M)
M= 0;



% ----------------------------- Eingabedaten -------------------------------

% Anfangstemperatur
% lexikographische Nummerierung der Gitterpunkte
% von x=0 bis x=xl
T0= zeros((mry+1)*(mrz+1), (mrx+1));

% Daten auf den beiden gegenueberliegenden
% Wuerfelseiten (x=0 & x=xl)
% lexikographische Nummerierung der Gitterpunkte
% fuer x=0 bzw. x=xl
S1= -4*ones((mry+1)*(mrz+1), ndt+1);
S2= 4*ones((mry+1)*(mrz+1), ndt+1);


% Loesung mit DROPS
sol= ipdrops(T0, S1, S2, M, xl, yl, zl, nu, mrx, mry, mrz, dt, theta, cgtol, cgiter);


function lambda=simipad(dTmess)

% gemeinsame Daten der direkten Probleme holen
[nu,xl,yl,zl,M,theta,cgtol,cgiter,Flag,BndRef,nis,ni,func]=getAdjointData;

npx= (ni(1)+1);
npyz= (ni(2)+1)*(ni(3)+1);
dt= 0.01;
ndt= 100;
time=[0:dt:ndt*dt];
%
%mr= 3;
%npc= 2^mr+1;
%nfp= (npc)^2;
%
T0= zeros(npyz, npx);
qh= dTmess;
qc= zeros(npyz, ndt+1);
%qc= qcp*ones(nfp, ndt+1);

% Loesung mit DROPS
cd ..;
T= ipdrops(T0, qh, qc, M, xl, yl, zl, nu, nis(1), nis(2), nis(3), dt, theta, cgtol, cgiter, Flag, BndRef);
cd optimization;
lambda=[T0(:,M+1) T];



%figure(2);plotsol(npc,ndt,sol)

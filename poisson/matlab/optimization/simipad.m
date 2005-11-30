function [Iter,lambda]=simipad(dTmess,C_data,L_data)

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
qh= dTmess;
qc= zeros(npyz, ndt+1);
%qc= qcp*ones(nfp, ndt+1);

% Loesung mit DROPS
work_dir= cd;
cd ..;
[Val1,T]= ipdrops(T0, qh, qc, M, xl, yl, zl, nu, nis(1), nis(2), nis(3), dt, theta, cgtol, cgiter, Flag, BndRef);
cd( work_dir);

Plane_Pos= round(M*nix/xl);
lambda= [T0(:,Plane_Pos+1) T];
Iter= Val1;

%figure(2);plotsol(npc,ndt,sol)

return

clear;

Nx=16; Ny=4; Nz=4;        % number of grid points
lx= 100; ly= 0.3; lz= 1;  % length of domain
Nt= 100; dt= 0.1;         % number and length of time steps

T_in= zeros(Ny*Nz, Nt);
T0  = zeros(Nx*Ny*Nz, 1);
F   = zeros(Nx*Ny*Nz, Nt);

Tsol= zeros(Nx*Ny*Nz, Nt-1);

[MaxIter, Tsol]= ipfilm(T0, T_in, F, 0.26, lx, ly, lz, Nx-1, Ny-1, Nz-1, dt, Nt-1, 0.5, 1e-10, 5000, 77);
  

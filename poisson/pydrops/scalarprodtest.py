import drops
import numpy as np

f = 1
Nx = 10*f
Ny = 10*f
Nz = 10*f
Nt = 1

lx = 1.
ly = 1.
lz = 1.
tmax = 0.
#dt = 0.1
#theta = 1.0
drops.setup_scalar_product_matrices(Nx,Ny,Nz,Nt,lx,ly,lz,tmax, False)

#print sum(retval[1:])

Xmesh, Ymesh, Zmesh, Tmesh = np.mgrid[0:Nx,0:Ny,0:Nz,0:Nt]
X = Xmesh.astype(float)/(Nx-1)*lx
Y = Ymesh.astype(float)/(Ny-1)*ly
#Z = Zmesh.astype(float)/(Nz-1)*lz
#T = Tmesh.astype(float)/(Nt-1)*tmax

#D = 100*ones([Nx, Ny, Nz, Nt])

trueval = 4/np.pi**2*lx*ly
c = np.sqrt(np.sin(X/lx*np.pi/2))*np.sqrt(np.sin((ly-Y)/ly*np.pi/2))

retval = drops.scalar_product(c,c)
err = np.abs(trueval - retval)
print "err", err[0]

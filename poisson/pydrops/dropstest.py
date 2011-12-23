import drops
from numpy import sin, pi, ones, zeros, mgrid
from matplotlib.pylab import *

f = 1
Nx = 10*f
Ny = 10*f
Nz = 10*f
Nt = 1

lx = 1.
ly = 1.
lz = 1.
tmax = 0.
dt = 0.1
theta = 1.0

drops.setup_scalar_product_matrices(Nx, Ny, Nz, Nt, lx, ly, lz, tmax, False)

Xmesh, Ymesh, Zmesh, Tmesh = mgrid[0:Nx,0:Ny,0:Nz,0:Nt]
X = Xmesh.astype(float)/(Nx-1)*lx
Y = Ymesh.astype(float)/(Ny-1)*ly

D = ones([Nx, Ny, Nz, Nt])

c = 1 + sin(X/lx*pi/2)*(1-sin((ly-Y)/ly*pi/2))

C0 = zeros([Nx, Ny, Nz])
b_in = c[0,:,:,:]
b_interface = c[:,-1,:,:]
source = -D *( -(pi/2/lx)**2*sin(pi/2*X/lx)*(1-sin(pi/2*(ly-Y)/ly)) + (pi/2/ly)**2*sin(pi/2*X/lx)*sin(pi/2*(ly-Y)/ly))

Dmol = 0.
uN = 0.
flag_pr = True
flag_bc = True
flag_supg = False

retval =  drops.convection_diffusion(C0,  b_in, source, D, b_interface, "dropstest.json")

dV  = lx*ly*lz/(Nx-1)/(Ny-1)/(Nz-1)
err = ((retval - c)*dV)**2
#print np.sqrt(np.sum(err))
if -1.2345 in c:
    print 'err'

diff = retval - c
print "err: ", diff
print "retval:", retval
sperr = "Scalar product: ", np.sqrt(drops.scalar_product(diff, diff))
print "vector norm: ", linalg.norm(err)
print sperr
#print retval
#print source

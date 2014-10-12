#! /usr/bin/python

# ...
try:
    from matplotlib import pyplot as plt
    PLOT=True
except ImportError:
    PLOT=False
# ...
import numpy                as np
from pigasus.fem.basicPDE import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

# ...
sin = np.sin ; pi = np.pi ; exp = np.exp
# ...

# ...
tc = {}

tc['A']  = lambda x,y : [1., 0., 0., 1.]

kx = 2.5 * pi ; ky = 1.5 * pi
xc = 4. / 5. ; yc = 1. / 4.

# ...
# exact solution
# ...
u = lambda x,y : [ sin ( kx * x - kx * xc ) * sin ( ky * y - ky * yc ) ]
#u = lambda x,y : [exp ( 0.5 * ( x**2 + y**2 ) )]
tc['u'] = u
# ...

# ...
# rhs
# ...
f = lambda x,y : [ ( kx**2 + ky**2 ) * sin ( kx * x - kx * xc ) * sin ( ky * y - ky * yc ) ]
#f = lambda x,y : [  -1.0*x**2*exp(0.5*x**2 + 0.5*y**2) - 1.0*y**2*exp(0.5*x**2 + 0.5*y**2) - 2.0*exp(0.5*x**2 + 0.5*y**2)]
tc['f'] = f
# ...

tc['Dirichlet'] = [[1]]
# ...
# values of u at the boundary
# ...
bc_dirichlet={}
bc_dirichlet [0,0] = u
#bc_dirichlet [0,1] = u
bc_dirichlet [0,2] = u
bc_dirichlet [0,3] = u

tc['bc_dirichlet'] = bc_dirichlet
# ...

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 31

try:
    ny = int(sys.argv[2])
except:
    ny = 31

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2
#-----------------------------------
# ...
from igakit.cad_geometry import square as domain
geo = domain(n=[nx,ny],p=[px,px])

# ...
PDE = basicPDE(geometry=geo, testcase=tc)
# ...

# ...
PDE.assembly()
PDE.solve()
# ...

# ...
normU = PDE.norm(exact=u)
print "norm U   = ", normU
# ...

# ...
if PLOT:
    PDE.plot()  ; plt.colorbar(); plt.title('$u_h$')
    plt.savefig(filename.split('.py')[0]+'.png', format='png')
    plt.clf()
# ...

PDE.free()

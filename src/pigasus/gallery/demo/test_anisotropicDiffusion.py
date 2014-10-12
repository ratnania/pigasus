#! /usr/bin/python

# ...
try:
    from matplotlib import pyplot as plt
    PLOT=True
except ImportError:
    PLOT=False
# ...
from pigasus.gallery.parabolic import onestep
from igakit.cad_geometry import square as domain
import numpy                as np
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

# ...
sin = np.sin ; pi = np.pi; exp = np.exp ; cos = np.cos
# ...

def tc(eps, dt, alpha, m1, m2, s1, s2):
    # ...
    testcase = {}

    scl = eps * dt * alpha
    phi = 2 * (1./3) * pi
    c = cos(phi)
    s = sin(phi)

    testcase['A'] = lambda x,y : [  scl * c**2 + s**2  \
                                 , ( scl - 1 ) * c * s \
                                 , ( scl - 1 ) * c * s \
                                 , scl * s**2 + c**2]

    testcase['b'] = lambda x,y: [1.]

    # ...
    def f(x,m,s):
        return exp(-(x-m)**2/(2.0*s))
    # ...

    # ...
    testcase['u'] = lambda x,y : [ 0. ]
    testcase['f'] = lambda x,y : [ f ( x , m1 , s1 ) * f ( y , m2 , s2 ) ]
    # ...

    return testcase

#-----------------------------------


#-----------------------------------
AllDirichlet = True

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

#-----------------------------------
geo = domain(n=[nx,ny],p=[px,py])

alpha   = 0.
dt      = 0.5
eps     = 1.e-3
niter   = 30
m1      = 0.5
m2      = 0.5
s1      = 0.01
s2      = 0.01

tc_E = tc(eps, dt,    alpha, m1, m2, s1, s2)
tc_I = tc(eps, dt, alpha-1., m1, m2, s1, s2)

# ...
PDE = onestep(geometry=geo, list_tc=[tc_E, tc_I], AllDirichlet=AllDirichlet)
# ...

# ...
PDE.assembly()
PDE.solve(niter)
# ...

# ...
if PLOT:
    PDE.plot() ; plt.colorbar(); plt.title('$u_h$')
    plt.savefig(filename.split('.py')[0]+'.png', format='png')
    plt.clf()
#...

PDE.free()

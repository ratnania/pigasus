#! /usr/bin/python

# ...
try:
    from matplotlib import pyplot as plt
    PLOT=True
except ImportError:
    PLOT=False
# ...
import numpy                as np
from pigasus.gallery.poisson import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

# ...
sin = np.sin ; cos = np.cos ; pi = np.pi ; exp = np.exp
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

from igakit.cad_geometry import quart_circle as domain
geo = domain(n=[nx,ny],p=[px,py])
#-----------------------------------

# ...
# exact solution
# ...
R = 1.
r = 0.5
c = 1. # for neumann
#c = pi / (R**2-r**2) # for all dirichlet bc
u = lambda x,y : [ x * y * sin ( c * (R**2 - x**2 - y**2 )) ]
# ...

# ...
# rhs
# ...
f = lambda x,y : [4*c**2*x**3*y*sin(c*(R**2 - x**2 - y**2)) \
                  + 4*c**2*x*y**3*sin(c*(R**2 - x**2 - y**2)) \
                  + 12*c*x*y*cos(c*(R**2 - x**2 - y**2)) ]
# ...

# ...
# values of gradu.n at the boundary
# ...
gradu   = lambda x,y : [-2*c*x**2*y*cos(c*(R**2 - x**2 - y**2)) + y*sin(c*(R**2
                                                                           -
                                                                           x**2
                                                                           -
                                                                           y**2)) \
                       ,-2*c*x*y**2*cos(c*(R**2 - x**2 - y**2)) + x*sin(c*(R**2 - x**2 - y**2)) ]

def func_g (x,y) :
    du  = gradu (x, y)
    return [ du[0] , du[1] ]
# ...

# ...
# values of u at the boundary
# ...

bc_neumann={}

bc_neumann [0,0] = func_g
Dirichlet = [[1,2,3]]

#AllDirichlet = True
# ...

# ...
try:
    bc_dirichlet
except NameError:
    bc_dirichlet = None
else:
    pass

try:
    bc_neumann
except NameError:
    bc_neumann = None
else:
    pass

try:
    AllDirichlet
except NameError:
    AllDirichlet = None
else:
    pass

try:
    Dirichlet
except NameError:
    Dirichlet = None
else:
    pass

try:
    Metric
except NameError:
    Metric = None
else:
    pass
# ...

# ...
PDE = poisson(geometry=geo, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann,
              AllDirichlet=AllDirichlet, Dirichlet=Dirichlet,metric=Metric)
# ...

# ...
PDE.assembly(f=f)
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

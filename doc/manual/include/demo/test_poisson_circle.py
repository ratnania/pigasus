#! /usr/bin/python

# ...
try:
    from matplotlib import pyplot as plt
    PLOT=True
except ImportError:
    PLOT=False
# ...
import numpy                as np
from igakit.cad_geometry import circle as domain
from igakit.cad_geometry import square as patch
from pigasus.gallery.poisson import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

# ...
sin = np.sin ; cos = np.cos
# ...

# ...
# exact solution
# ...
u = lambda x,y : [sin ( 1.0 - x**2 - y**2 ) ]
# ...

# ...
# rhs
# ...
f = lambda x,y : [4.0 * ( x**2 + y**2 ) * sin ( 1.0 - x**2 - y**2 ) + 4.0 * cos ( 1.0 - x**2 - y**2 ) ]
# ...

#-----------------------------------

#-----------------------------------
# ...
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


geo   = domain(n=[nx,ny],p=[px,py])
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

PDE = poisson(geometry=geo, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann,
              AllDirichlet=AllDirichlet, Dirichlet=Dirichlet,metric=Metric)

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

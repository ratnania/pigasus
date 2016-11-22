# -*- coding: UTF-8 -*-
#! /usr/bin/python
from pigasus.utils.manager import context

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
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ...
sin = np.sin ; cos = np.cos ; pi = np.pi ; exp = np.exp ; sqrt = np.sqrt
# ...

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 1

try:
    ny = int(sys.argv[2])
except:
    ny = 1

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2

from caid.cad_geometry import circle as domain
radius = 0.5
geo = domain(radius=radius,n=[nx,ny],p=[px,px])
#-----------------------------------

# ...
# exact solution
# ...
meanU = 0.60125862034
u = lambda x,y : [sin ( 1.0 - x**2 - y**2 )-meanU/(pi*radius**2)]
# ...

# ...
# rhs
# ...
f = lambda x,y : [4.0 * ( x**2 + y**2 ) * sin ( 1.0 - x**2 - y**2 ) + 4.0 * cos ( 1.0 - x**2 - y**2 ) ]
# ...

# ...
# values of gradu.n at the boundary
# ...
gradu   = lambda x,y : [ -2. * x * cos ( 1.0 - x**2 - y**2 ) \
                       , -2. * y * cos ( 1.0 - x**2 - y**2 ) ]

def func_g (x,y) :
    du  = gradu (x, y)
    return [ du[0] , du[1] ]
# ...

# ...
# values of u at the boundary
# ...
bc_neumann={}
bc_neumann [0,0] = func_g
bc_neumann [0,1] = func_g
bc_neumann [0,2] = func_g
bc_neumann [0,3] = func_g
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

with context():
    # ...
    PDE = poisson(geometry=geo, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann,
                  AllDirichlet=AllDirichlet, Dirichlet=Dirichlet,metric=Metric)
    # ...

    # ...
    PDE.assembly(f=f)
    PDE.solve()
    from scipy.io import mmwrite
    mmwrite("M.mtx", PDE.system.get())
    np.savetxt("rhs.txt", PDE.rhs.get())
    # ...

    # ...
    normU = PDE.norm(exact=u)
    print("norm U   = ", normU)
    # ...

    # ...
    if PLOT:
        PDE.plot()  ; plt.colorbar(); plt.title('$u_h$')
        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()
    # ...

    PDE.free()

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
from caid.cad_geometry import square as domain
from pigasus.gallery.poisson import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ...
sin = np.sin ; pi = np.pi
# ...

# ...
kx = 2. * pi
ky = 2. * pi

# ... exact solution
u = lambda x,y : [sin ( kx * x ) * sin ( ky * y )]
# ... rhs
f = lambda x,y : [( kx**2 + ky**2 ) * sin ( kx * x ) * sin ( ky * y )]


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


geo = domain(n=[nx,ny],p=[px,py])
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
    PDE = poisson(geometry=geo, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann,
                  AllDirichlet=AllDirichlet, Dirichlet=Dirichlet,metric=Metric)

    # ...

    PDE.assembly(f=f)
    M = PDE.system
    M.save("M.mtx")
    PDE.solve()
    # ...

    # ...
    normU = PDE.norm(exact=u)
    print("norm U   = ", normU)
    # ...

    U = PDE.unknown
    U.export("U.pfl")

    # ...
    if PLOT:
        PDE.plot()  ; plt.colorbar(); plt.title('$u_h$')
        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()
    # ...

    PDE.free()

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
from caid.cad_geometry import line as domain
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

# ... exact solution
u = lambda x : [sin ( kx * x )]
# ... rhs
f = lambda x : [( kx**2) * sin ( kx * x )]

#-----------------------------------
nx      = 63
px      = 2
AllDirichlet = True

geo = domain(n=[nx],p=[px])
#-----------------------------------


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
    # ...

    # ...
    normU = PDE.norm(exact=u)
    print("norm U   = ", normU)
    # ...

    # ...
    if PLOT:
        PDE.plot()  ; plt.title('$u_h$')
        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()
    # ...

    PDE.free()

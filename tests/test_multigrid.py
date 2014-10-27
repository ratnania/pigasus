# -*- coding: UTF-8 -*-
#! /usr/bin/python
from pigasus.utils.manager import context

from pigasus.gallery.multigrid import *
from pigasus.gallery.poisson import *
from pigasus.utils.utils import hierarchical_geometries
from caid.cad_geometry import square as domain
import numpy                as np
import time
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

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

try:
    nlevel = int(sys.argv[5])
except:
    nlevel = 4
#-----------------------------------

geo = domain(n=[nx,ny],p=[px,py])
list_geometry = hierarchical_geometries(geo, nlevel, domain)
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
geo_h = list_geometry[-1]

with context():

    PDE = poisson(geometry=geo_h, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann \
                  , AllDirichlet=AllDirichlet, Dirichlet=Dirichlet,metric=Metric)
    rhs = PDE.rhs
    # ...

    # ...
    MG = multigrid(PDE, list_geometry=list_geometry)
    PDE.assembly()
    # ...

    # ...
    mg_residuals = MG.solve(rhs, verbose=True, accel=None)
    mg_residuals = MG.solve(rhs, verbose=True, accel='gmres')
    # ...

    PDE.free()

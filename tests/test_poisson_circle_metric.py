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
    nx = 3

try:
    ny = int(sys.argv[2])
except:
    ny = 3

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2

from caid.cad_geometry import circle as domain
from caid.cad_geometry import cad_geometry, cad_nurbs
from pigasus.fem.metric import metric
from caid.cad_geometry import square as patch
from caid.cad_geometry import square

TYPE = None
#TYPE = "mapping"
#TYPE = "points"
#TYPE = "analytic"

if TYPE is None:
    geo   = domain(n=[nx,ny],p=[px,py])

# -----------------
if TYPE is not None:
#    geo = cad_geometry("domain.xml")

#    geo   = patch(n=[nx,ny],p=[px,py])

    geo_s = square(p=[2,2])
    nrb = geo_s[0]
    U,V = nrb.knots
    #C   = nrb.points
    C   = np.zeros_like(nrb.points)

    s = 1./np.sqrt(2)
    weights         = np.ones((3,3))
    weights[1,0]    = s
    weights[0,1]    = s
    weights[2,1]    = s
    weights[1,2]    = s

    srf = cad_nurbs([U,V], C, weights=weights)
    geo = cad_geometry()
    geo.append(srf)
    geo._internal_faces = geo_s._internal_faces
    geo._external_faces = geo_s._external_faces
    geo._connectivity   = geo_s._connectivity
    #
    srf = geo[0]
    dim = srf.dim
    n = [nx,ny]
    list_t = []
    for axis in range(0,dim):
        ub = srf.knots[axis][0]
        ue = srf.knots[axis][-1]
        list_t.append(np.linspace(ub,ue,n[axis]+2)[1:-1])

    p =[px,py]
    list_p = []
    for axis in range(0,dim):
        list_p.append(p[axis] - srf.degree[axis])

    geo.refine(list_t=list_t, list_p=list_p)
# -----------------


if TYPE is None:
    Metric = None
    geo   = domain(n=[nx,ny],p=[px,py])

if TYPE == "mapping":
    geo_m = domain()
    Metric = metric(geometry=geo_m, with_igakit=False)

if TYPE == "analytic":
    F  = lambda r,t : [r * cos(2. * pi * t), r * sin(2. * pi * t)]
    DF = lambda r, t : [cos(2. * pi * t) \
                       , - 2. * pi * r * sin(2. * pi * t) \
                       , sin(2. * pi * t) \
                       , 2. * pi * r * cos(2. * pi * t)]
    Metric = metric(analytic=[F,DF])

if TYPE == "points":
    lpi_shape = np.asarray([int(x) for x in np.genfromtxt('shape.txt')])
#    lpr_points = np.genfromtxt('points_adv.txt').reshape(lpi_shape)
    lpr_points = np.genfromtxt('pts_adv.txt').reshape(lpi_shape)
    Metric = metric(points = lpr_points)
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

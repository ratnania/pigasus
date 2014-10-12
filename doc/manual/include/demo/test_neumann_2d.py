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

TEST = 0

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

#from igakit.cad_geometry import circle as domain
from igakit.cad_geometry import square as domain
geo = domain(n=[nx,ny],p=[px,px])
#execfile('gen_square_def.py')
#-----------------------------------

# ...
if TEST==0:
    kx = 1./3. * pi ; ky = 5./3. * pi
    xc = 1./3. ; yc = 1./4.
    meanU = 0.030590513542813665

if TEST==1:
    kx = pi ; ky =  pi
    xc = 0.5  ; yc = 0.5
    meanU = 0.0

if TEST==2:
    kx = pi ; ky =  0.25 * pi
    xc = 0.  ; yc = 0.

if TEST==21:
    kx = 1.25 * pi ; ky =  pi
    xc = 0.  ; yc = 0.

if TEST==3:
    kx = 1.25 * pi ; ky =  0.25 * pi
    xc = 0.  ; yc = 0.

if TEST==4:
    kx = 1.25 * pi ; ky =  0.25 * pi
    xc = 0.  ; yc = 0.3

# ...
# exact solution
# ...
u = lambda x,y : [ sin ( kx * x - kx * xc ) * sin ( ky * y - ky * yc ) ]
if TEST in [0,1]:
    u = lambda x,y : [ sin ( kx * x - kx * xc ) * sin ( ky * y - ky * yc ) - meanU]
# ...

# ...
# rhs
# ...
f = lambda x,y : [ ( kx**2 + ky**2 ) * sin ( kx * x - kx * xc ) * sin ( ky * y - ky * yc ) ]
# ...

# ...
# values of gradu.n at the boundary
# ...
gradu   = lambda x,y : [kx*sin(ky*y - ky*yc)*cos(kx*x - kx*xc) \
                        , ky*sin(kx*x - kx*xc)*cos(ky*y - ky*yc)]

def func_g (x,y) :
    du  = gradu (x, y)
    return [ du[0] , du[1] ]
# ...

# ...
# values of u at the boundary
# ...
if TEST in [0,1]:
    bc_neumann={}
    bc_neumann [0,0] = func_g
    bc_neumann [0,1] = func_g
    bc_neumann [0,2] = func_g
    bc_neumann [0,3] = func_g

if TEST == 2:
    bc_neumann={}
    bc_neumann [0,2] = func_g

    Dirichlet = [[0,1,3]]


if TEST == 21:
    bc_neumann={}
    bc_neumann [0,3] = func_g

    Dirichlet = [[0,1,2]]

if TEST == 3:
    bc_neumann={}
    bc_neumann [0,2] = func_g
    bc_neumann [0,3] = func_g

    Dirichlet = [[0,1]]

if TEST == 4:
    bc_neumann={}
    bc_neumann [0,0] = func_g
    bc_neumann [0,2] = func_g
    bc_neumann [0,3] = func_g

    Dirichlet = [[1]]
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

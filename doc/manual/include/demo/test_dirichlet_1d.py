#! /usr/bin/python

# ...
try:
    from matplotlib import pyplot as plt
    PLOT=True
except ImportError:
    PLOT=False
# ...
import numpy as np
from pigasus.fem.basicPDE import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

# ...
sin = np.sin ; pi = np.pi
# ...

# ...
tc = {}

tc['A']  = lambda x : [1.]

kx = 0.5 * pi

# ...
# exact solution
# ...
u = lambda x : [ sin ( kx * x ) ]
tc['u'] = u
# ...

# ...
# rhs
# ...
f = lambda x : [ ( kx**2) * sin ( kx * x ) ]
tc['f'] = f
# ...

tc['Dirichlet'] = [[0]]
# ...
# values of u at the boundary
# ...
bc_dirichlet={}
bc_dirichlet [0,1] = u

tc['bc_dirichlet'] = bc_dirichlet
# ...

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 31

try:
    px = int(sys.argv[2])
except:
    px = 2
#-----------------------------------

# ...
from igakit.cad_geometry import line as domain
geo = domain(n=[nx],p=[px])
# ...

# ...
PDE = basicPDE(geometry=geo, testcase=tc)
# ...

# ...
PDE.assembly()
PDE.solve()

# ...
normU = PDE.norm(exact=u)
print "norm U   = ", normU
# ...

# ...
if PLOT:
    PDE.plot()  ; plt.title('$u_h$')
    plt.savefig(filename.split('.py')[0]+'.png', format='png')
    plt.clf()
# ...

PDE.free()

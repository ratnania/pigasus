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
import numpy as np
from pigasus.fem.basicPDE import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ...
sin = np.sin ; pi = np.pi
# ...

# ...
tc = {}

tc['A']  = lambda x : [1.]

kx = 2. * pi

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

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 5

try:
    px = int(sys.argv[2])
except:
    px = 3
#-----------------------------------

# ...
from caid.cad_geometry import line as domain
geo = domain(n=[nx],p=[px], periodic=True)
# ...

with context():
    # ...
    PDE = basicPDE(geometry=geo, testcase=tc)
    # ...

    # ...
    PDE.assembly()
    M = PDE.system
    M.save("M.mtx")
    PDE.solve()

    # ...
    normU = PDE.norm(exact=u)
    print("norm U   = ", normU)
    # ...

    # ...
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown
    U.export(filename.split('.py')[0]+".pfl")
    # ...

    # ...
#    if PLOT:
#        PDE.plot()  ; plt.title('$u_h$')
#        plt.savefig(filename.split('.py')[0]+'.png', format='png')
#        plt.clf()
    # ...

    PDE.free()

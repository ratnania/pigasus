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
from pigasus.gallery.bilaplacian import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ...
sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt ; pi = np.pi
# ...

#-----------------------------------
kx = 2. * pi ; ky = 2. * pi

# exact solution
u = lambda x,y : [sin ( kx * x ) * sin ( ky * y )]

# rhs
f = lambda x,y : [( kx**4 + ky**4 ) * sin ( kx * x ) * sin ( ky * y )]
#-----------------------------------

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
    px = 3

try:
    py = int(sys.argv[4])
except:
    py = 3
#-----------------------------------

#-----------------------------------
# ...
geo = domain(n=[nx,ny], p=[px,py])

with context():
    PDE = bilaplacian(geometry=geo)
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

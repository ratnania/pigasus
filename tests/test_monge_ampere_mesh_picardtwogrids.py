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
from caid.cad_geometry import square
from caid.cad_geometry import circle
from caid.cad_geometry import quart_circle
from caid.cad_geometry import annulus
import numpy                as np
from time import time
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ... import picard from monge_ampere module
from pigasus.utils.load import load
monge_ampere    = load("monge_ampere")
picardTwoGrids  = monge_ampere.picardTwoGrids
testcase        = monge_ampere.testcase
# ...

abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 15

try:
    ny = int(sys.argv[2])
except:
    ny = 15

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2

geo   = square(n=[nx,ny], p=[px,py])
#geo   = circle(radius=1.,n=[nx,ny], p=[px,py])
#geo   = quart_circle(n=[nx,ny], p=[px,py])
#geo   = annulus(n=[nx,ny], p=[px,py])

#from caid.cad_geometry import cad_geometry as domain
#geo = domain("input/iter_inner.xml")
#-----------------------------------

#-----------------------------------
tc = testcase(1)
rho0 = tc.rho0
rho1 = tc.rho1

verbose = 1

#    p_H = [ 5, 5]
#    p_h = [ 5, 5]

p_H = [ 3, 3]
p_h = [ 3, 3]

#    p_H = [ 2, 2]
#    p_h = [ 2, 2]

# TEST 1
#    # p = 5
#    rtol_H = 1.e-4
#    rtol_h = 1.e-8

# p = 3
#    rtol_H = 1.e-4
#    rtol_h = 1.e-6

# p = 2
#    rtol_H = 1.e-4
#    rtol_h = 1.e-4

# TEST 3
# p = 3, 5
rtol_H = 1.e-3
rtol_h = 1.e-3

#    # p = 2
##    rtol_H = 1.e-2
##    rtol_h = 1.e-2


rtol2_H = 1.e-6
#    rtol2_h = 1.e-6
rtol2_h = 1.e-9

maxiter_H = 40
maxiter_h = 40

n_H = [7,7]
nstage =  1
#-----------------------------------

with context():

    # ...
    n_h = []
    for axis in range(0,2):
        n = n_H[axis]
        for i in range(0, nstage):
            n = 2*n+1
        n_h.append(n)

    print(">>>> coarse grid ", n_H, " with splines of degree ", p_H)
    print(">>>> fine   grid ", n_h, " with splines of degree ", p_h)

    geo_H = square(n=n_H, p=p_H)
    geo_h = square(n=n_h, p=p_h)
    # ...

    # ...
    # values of gradu.n at the boundary
    # ...
    def func_g(x,y):
        return [x,y]
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
    PDE_h = picardTwoGrids(geometry=geo_h, geometry_H=geo_H, bc_neumann=bc_neumann)
    # ...

    # ...
    Errors_h, ErrorsH1_h, Errors_H, ErrorsH1_H = \
            PDE_h.solve(  rho0, rho1, c_rho=None, u0=None \
                        , maxiter=[maxiter_H,maxiter_h] \
                        , rtol=[rtol_H,rtol_h] \
                        , rtol2=[rtol2_H,rtol2_h] \
                        , verbose=1)
    # ...

    # ...
    if PLOT:
        PDE_h.plotMesh(ntx=60, nty=60)
        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()
    # ...

    np.savetxt("Errors.txt", np.asarray(Errors_h))

    PDE_h.free()

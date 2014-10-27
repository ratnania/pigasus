#! /usr/bin/python

import sys
import numpy as np
from pigasus.gallery.basicPDE import *

import caid.cad_geometry  as cg
import matplotlib.pyplot    as pl
import numpy                as np

# ...
sin = np.sin ; cos = np.cos ; pi = np.pi
# ...

# ...
def func_n (x, y):
    list_nx = []
    list_ny = []
    for (u,v) in zip(x,y):
        nx = 0. ; ny = 0.
        if np.allclose(u, 0.):
            nx = -1.
            ny = 0.
        if np.allclose(u, 1.):
            nx = 1.
            ny = 0.

        if np.allclose(v, 0.):
            nx = 0.
            ny = -1.
        if np.allclose(v, 1.):
            nx = 0.
            ny = 1.

        list_nx.append(nx)
        list_ny.append(ny)

    nx = np.asarray(list_nx)
    ny = np.asarray(list_ny)

    return nx, ny
# ...

def testcase():
    # ...
    testcase = {}
    testcase['Dirichlet'] = [[0,1]]

    testcase['A']  = lambda x,y : [1., 0., 0., 1.]

    xc = 0.
    yc = 0.

    kx = .5 * pi
    ky = .5 * pi

    # ...
    # exact solution
    # ...
    u = lambda x,y : [ sin(kx*(x-xc)) * sin(ky*(y-yc)) ]
    testcase['u'] = u
    # ...

    # ...
    # rhs
    # ...
    f = lambda x,y : [ (kx**2 + ky**2) * sin(kx*(x-xc)) * sin(ky*(y-yc)) ]
    testcase['f'] = f
    # ...

    # ...
    # values of gradu.n at the boundary
    # ...
    gradu   = lambda x,y : [  kx * cos(kx*(x-xc)) * sin(ky*(y-yc)) \
                                , ky * sin(kx*(x-xc)) * cos(ky*(y-yc)) ]

    def g (x,y) :
        du  = gradu (x, y)
        nx, ny = func_n (x, y)

        return [nx * du[0] + ny * du[1]]

    testcase['g'] = g
    # ...

    # ...
    # values of u at the boundary
    # ...
    bc_neumann={}
#    bc_neumann [0,0] = g
#    bc_neumann [0,1] = g
    bc_neumann [0,2] = g
#    bc_neumann [0,3] = g
#
    testcase['bc_neumann'] = bc_neumann
    # ...

    # ...
    # values of u at the boundary
    # ...
    bc_dirichlet={}
    bc_dirichlet [0,3] = u

    testcase['bc_dirichlet'] = bc_dirichlet
    # ...

    return testcase
#-----------------------------------

#-----------------------------------
# ...
from caid.cad_geometry import square as domain
tc = testcase()
geo = domain(n=[31,31],p=[2,2])
PDE = basicPDE(geometry=geo, testcase=tc)
# ...

# ...
PDE.assembly()
PDE.solve()
normU = PDE.norm()
print "norm U = ", normU
PDE.plot()  ; pl.show()
# ...

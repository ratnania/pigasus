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
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

from pigasus.fem.block_basicPDE import *
from caid.cad_geometry import square as domain
from numpy import pi, sin

with context():
    nx = 15 ; ny = 15
    px = 2 ; py = 2
    geo = domain(n=[nx,ny],p=[px,px])

    #-----------------------------------
    def testcase():
        # implicit part
        tc = {}

        kx = 2. * pi
        ky = 2. * pi

        # ... exact solution
        u = lambda x,y : [sin ( kx * x ) * sin ( ky * y )]
        # ... rhs
        f = lambda x,y : [( kx**2 + ky**2 ) * sin ( kx * x ) * sin ( ky * y )]

        A = lambda x,y : [  1., 0. \
                          , 0., 1. ]

        tc['u']  = u
        tc['f']  = f
        tc['A']  = A

        tc['AllDirichlet'] = True

        return tc
    #-----------------------------------

    size = [2,2]

    dict_testcases = {}
    dict_testcases[0,0] = testcase()
#    dict_testcases[0,1] = testcase()
#    dict_testcases[1,0] = testcase()
    dict_testcases[1,1] = testcase()

    PDEs = block_basicPDE(size, dict_testcases=dict_testcases, geometry=geo)
    PDEs.assembly()

#    rhs = PDEs.rhs.get()
#    system = PDEs.system
#    matrix = system.get()
#
#    Y = matrix.dot(rhs)
#    print Y.shape

    rhs = PDEs.rhs
    PDEs.solve(rhs)

    # ... example of multiplication with a list of fields
    unknowns = PDEs.unknowns
    system = PDEs.system
    X = system.dot(unknowns)
    print([np.linalg.norm(x.get()-r.get()) for (x,r) in zip(X,rhs)])
    # ...

    # ... example of multiplication with a list of numpy arrays
    unknowns = [U.get() for U in PDEs.unknowns]
    system = PDEs.system
    X = system.dot(unknowns)
    print([np.linalg.norm(x-r.get()) for (x,r) in zip(X,rhs)])
    # ...


#    for U in PDEs.unknowns:
#        print U.get()

    print(PDEs.norms())


    PDEs.free()

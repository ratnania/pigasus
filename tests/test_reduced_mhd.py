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
from pigasus.fem.basicPDE import *
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
#sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ...
sin = np.sin ; pi = np.pi ; exp = np.exp
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
#-----------------------------------

#-----------------------------------
def testcase_1():
    # implicit part
    tc = {}

    tc['b']  = lambda x,y : [1.]

    tc['AllDirichlet'] = True

    tc_i = tc

    # explicit part
    tc = {}

    tc['b']  = lambda x,y : [1.]
    tc['v']  = lambda x,y : [0., 0.]

    tc['AllDirichlet'] = True

    tc_e = tc

    return tc_i, tc_e
#-----------------------------------

#-----------------------------------
def testcase_2():
    # implicit part
    tc = {}

    tc['b']  = lambda x,y : [1.]

    tc['AllDirichlet'] = True

    tc_i = tc

    # explicit part
    tc = {}

    tc['b']  = lambda x,y : [1.]
    tc['v']  = lambda x,y : [0., 0.]

    tc['AllDirichlet'] = True

    tc_e = tc

    return tc_i, tc_e
#-----------------------------------

#-----------------------------------

#-----------------------------------
def func_field(U,x,y):
    # ...
    D    = U.evaluate(nderiv=1, parametric=False)

    _U   = D[0,0,:]
    Udx  = D[0,1,:]
    Udy  = D[0,2,:]

    return [Udy, -Udx]
# ...
#-----------------------------------

#-----------------------------------
from pigasus.fem.utils import function
def create_function(U):
    return function(func_field, fields=[U])
#-----------------------------------

#-----------------------------------
def get_unknwon(PDE):
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown

    return U
#-----------------------------------

# ...
tc_1i, tc_1e = testcase_1()
tc_2i, tc_2e = testcase_2()
# ...

# ...
from caid.cad_geometry import square as domain
geo = domain(n=[nx,ny],p=[px,px])

with context():

    # ...
    PDE_1e = basicPDE(geometry=geo, testcase=tc_1e)
    PDE_1i = basicPDE(geometry=geo, testcase=tc_1i, V=PDE_1e.V)
    PDE_2e = basicPDE(geometry=geo, testcase=tc_2e, V=PDE_1e.V)
    PDE_2i = basicPDE(geometry=geo, testcase=tc_2i, V=PDE_1e.V)
    # ...

    # ...
    U_1e = get_unknwon(PDE_1e)
    U_1i = get_unknwon(PDE_1i)
    U_2e = get_unknwon(PDE_2e)
    U_2i = get_unknwon(PDE_2i)
    # ...

    # ...
    Advection_1e = PDE_1e.advection
    Advection_1e.func = create_function(U_2e)
    # ...

    # ...
    Advection_2e = PDE_2e.advection
    Advection_2e.func = create_function(U_1e)
    # ...

    # ...
    PDE_1e.assembly()
    PDE_1i.assembly()
    PDE_2e.assembly()
    PDE_2i.assembly()
    # ...



    # ...
#    PDE_1.solve()
#    PDE_2.solve()
    # ...

    # ...
#    normU = PDE.norm(exact=u)
#    print "norm U   = ", normU
    # ...

    # ...
    PDE_1e.free()
    PDE_1i.free()
    PDE_2e.free()
    PDE_2i.free()
    # ...

# ...


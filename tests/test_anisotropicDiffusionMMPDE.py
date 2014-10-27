#! /usr/bin/python

import sys
import numpy as np

import caid.cad_geometry  as cg
import matplotlib.pyplot    as plt
import numpy                as np
from __main__ import __file__ as filename

# ...
sin = np.sin ; pi = np.pi; exp = np.exp ; cos = np.cos
# ...

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
#-----------------------------------

#-----------------------------------
def testcase(eps, dt, alpha, m1, m2, s1, s2):
    # ...
    testcase = {}
    testcase['AllDirichlet'] = True
#    testcase['Dirichlet'] = [[0,1,2,3]]

    scl = eps * dt * alpha
    phi = 2 * (1./3) * pi
    c = cos(phi)
    s = sin(phi)

    testcase['A'] = lambda x,y : [  scl * c**2 + s**2  \
                                 , ( scl - 1 ) * c * s \
                                 , ( scl - 1 ) * c * s \
                                 , scl * s**2 + c**2]

    testcase['b'] = lambda x,y: [1.]

    testcase['v'] = lambda x,y: [0., 0.] # will be defined after

    # ...
    def f(x,m,s):
        return exp(-(x-m)**2/(2.0*s))
    # ...

    # ...
    testcase['u'] = lambda x,y : [ 0. ]
    testcase['f'] = lambda x,y : [ f ( x , m1 , s1 ) * f ( y , m2 , s2 ) ]
    # ...

    return testcase
#-----------------------------------

#-----------------------------------
alpha   = 0.
beta    = 0.
t       = 0.
dt      = 0.5
eps     = 1.e-3
niter   = 30
m1      = 0.5
m2      = 0.5
s1      = 0.01
s2      = 0.01

tc_E = testcase(eps, dt,    alpha, m1, m2, s1, s2)
tc_I = testcase(eps, dt, alpha-1., m1, m2, s1, s2)
#-----------------------------------

# ...
from caid.cad_geometry import square as domain
geo = domain(n=[nx,ny],p=[px,px])
execfile('../run_anisotropicDiffusionMMPDE.py')
# ...

# ...
plt.savefig(filename.split('.py')[0]+'.png', format='png')
# ...

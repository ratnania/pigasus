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
from caid.cad_geometry import circle
from caid.cad_geometry import circle_5mp
from pigasus.gallery.poisson_nonlin import poisson_newton
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

exp = np.exp ; log = np.log ; sqrt = np.sqrt

#-----------------------------------
AllDirichlet = True

try:
    nx = int(sys.argv[1])
except:
    nx = 7

try:
    ny = int(sys.argv[2])
except:
    ny = 7

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2

radius = 1. / sqrt (2)
#geo = circle (radius=radius, n=[nx, ny], p=[px, py])
rmin = 0.5 * radius ; rmax = radius
geo = circle_5mp (rmin=rmin, rmax=rmax, n =[nx, ny], p =[px, py])
#-----------------------------------

# ...
u_exact = lambda x,y : [- 2.0 * log ( x**2 + y**2 + 0.5 )]

def F(U,x,y):
    _U = U.evaluate()
    return [4. * exp (_U)]

def dF (U,x, y):
    U = U.evaluate()
    return[-4 * exp (U)]
# ...

with context():

    PDE_newton = poisson_newton(  geometry=geo \
                         , AllDirichlet=AllDirichlet )

    print(">>> Solving using Newton <<<")
    # ...
    PDE = PDE_newton
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown
    # ...

    # ...
    list_L2, list_H1 = PDE_newton.solve(F, dF, u0=None, maxiter=40, rtol=1.e-6, verbose=True)

    print("norm using Newton  ", PDE_newton.norm(exact=u_exact))

    # ...
    if PLOT:
        fig = plt.figure()

        plt.subplot(121, aspect='equal')
        U.plot(withpcolor=True, n=[50,50]) ; plt.colorbar(orientation='horizontal') ; plt.title('$u_h$')


        # plot error evolution
        plt.subplot(122)
        plt.plot(list_L2, '-vb', label='$L^2$ norm')
        plt.plot(list_H1, '-xr', label='$H^1$ norm')
        plt.xlabel('N')
        plt.semilogy()
        plt.title('Norm evolution of $u^{n+1} - u^n$')
        plt.legend()
        # ...

        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()
    # ...

    PDE.free()

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
from pigasus.gallery.poisson_nonlin import poisson_picard
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

geo = circle (radius = 1. / sqrt (2), n =[nx, ny], p =[px, py])
#-----------------------------------

# ...
u_exact = lambda x,y : [- 2.0 * log ( x**2 + y**2 + 0.5 )]
# ...

with context():

    PDE = poisson_picard(  geometry=geo \
                         , AllDirichlet=AllDirichlet )

    # ...
    print(">>> Solving using Picard <<<")
    # ...
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown
    # ...

    from pigasus.fem.utils import function
    # ...
    def func(U,x,y):
        _U = U.evaluate()
        return [4. * exp (_U)]
    # ...

    F = function(func, fields=[U])

    list_L2, list_H1 = PDE.solve(F, u0=None, maxiter=50, rtol=1.e-6,
                                        verbose=True)

    print("norm using Picard  ", PDE.norm(exact=u_exact))

    # ...
    if PLOT:
        fig = plt.figure()

        plt.subplot(121, aspect='equal')
        U.fast_plot() ; plt.colorbar(orientation='horizontal') ; plt.title('$u_h$')

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

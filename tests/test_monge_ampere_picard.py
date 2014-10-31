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
from pigasus.gallery.poisson_nonlin import poisson_picard
from caid.cad_geometry import square
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

from time import time

exp = np.exp ; log = np.log ; sqrt = np.sqrt

#-----------------------------------
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
    px = 3

try:
    py = int(sys.argv[4])
except:
    py = 3
#-----------------------------------

#-----------------------------------
# ...
# exact solution
# ...
u_exact = lambda x,y : [exp ( 0.5 * ( x**2 + y**2 ) )]
u0      = lambda x,y : 0.5 * ( x**2 + y**2 )
f0      = lambda x,y : ( 1. + x**2 + y**2 ) * exp ( x**2 + y**2 )
# ...
#-----------------------------------

#-----------------------------------
geo = square (n =[nx, ny], p =[px, py])
#-----------------------------------

#-----------------------------------
# ...
# values of u at the boundary
# ...
bc_dirichlet={}
for data in geo.external_faces:
    patch_id = int(data[0]) ; face_id = int(data[1])
    bc_dirichlet[patch_id,face_id] = u_exact
# ...
#-----------------------------------

with context():
    PDE = poisson_picard(  geometry=geo \
                         , bc_dirichlet=bc_dirichlet)

    # ...
    print(">>> Solving using Picard <<<")
    # ...
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown

    V = PDE.space
    V.nderiv_pts = 2
    PDE.W.nderiv_pts = 2
    # ...

    def func(U,x,y):
        # ...
        D    = U.evaluate(nderiv=2, parametric=False)

        _U   = D[0,0,:]
        Udx  = D[0,1,:]
        Udy  = D[0,2,:]
        Udxx = D[0,3,:]
        Udxy = D[0,4,:]
        Udyy = D[0,5,:]

        f_values = f0(x,y)
        _F = - np.sqrt ( Udxx**2 + Udyy**2 + 2 * Udxy**2 + 2 * f_values )

        return [_F]
    # ...
    from pigasus.fem.utils import function
    F = function(func, fields=[U])

    list_L2, list_H1 = PDE.solve(F, u0=None, maxiter=300, rtol=1.e-6, verbose=True)

    u_exact = lambda x,y : [exp ( 0.5 * ( x**2 + y**2 ) )]
    print("norm using Picard  ", PDE.norm(exact=u_exact))

    # ...
    if PLOT:
        fig = plt.figure()

    #    plt.subplot(121, aspect='equal')
    #    U.fast_plot() ; plt.colorbar(orientation='horizontal') ; plt.title('$u_h$')
    #
    #    # plot error evolution
    #    plt.subplot(122)
    #    plt.plot(list_L2, '-vb', label='$L^2$ norm')
    #    plt.plot(list_H1, '-xr', label='$H^1$ norm')
    #    plt.xlabel('N')
    #    plt.semilogy()
    #    plt.title('Norm evolution of $u^{n+1} - u^n$')
    #    plt.legend()
    #    # ...
    #
    #    plt.savefig(filename.split('.py')[0]+'.png', format='png')
    #    plt.clf()
    # ...

        U = PDE.unknown_dirichlet
        tx = np.linspace(0.,1.,100)
        ty = np.linspace(0.,1.,200)
        uh = U(tx,ty, patch_id=0)
        nrb = geo[0]
        P = nrb(u=tx,v=ty)
        x = P[:,:,0]
        y = P[:,:,1]
    #    u = exp ( 0.5 * ( x**2 + y**2 ) )
        u = 0.
        plt.contourf(x,y,u-uh)  ; plt.colorbar(); plt.title('$u-u_h$')
        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()

    PDE.free()

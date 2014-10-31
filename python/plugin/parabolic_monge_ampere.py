# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from pigasus.gallery.basicPDE import *
import matplotlib.pyplot    as plt
import numpy                as np
from caid.cad_geometry import cad_nurbs
from __main__ import __file__ as filename

# ...
abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; log = np.log ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)
# ...

#-----------------------------------
n = [15,15]
p = [ 3, 3]
#-----------------------------------

#-----------------------------------
# ...
C0   = 1.0
rho0 = lambda x,y : 1.

#    C1   = 0.616805883732
#    t = 0.5
#    rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-t)**2+(y-0.5)**2-0.09)))

C1   =  1.75484181939
rho1 = lambda x,y : ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))

# ... test7
#xc = 0.7 ; yc = 0.5
#C1 = 0.281648379406
#
#r = lambda s,t : sqrt( (s-xc)**2 + (t-yc)**2 )
#theta = lambda s,t : atan(t-yc,s-xc)
#def rho1(s,t):
#    r_ = r(s,t) ;  t_ = theta(s,t)
#    val = C1 * (1. + 9./(1. + (10*r_*cos(t_-20*r_**2))**2) )
#    return val
# ...

c_rho = C0/C1
# ...
#-----------------------------------

#-----------------------------------
# ...
from caid.cad_geometry import square as domain
geo = domain(n=n,p=p)

gamma   = 11.
eps     = 1.
dt      = 1.
rtol    = 1.e-3
maxiter = 10
verbose = True

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
tc = {}
tc['A'] = lambda x,y : [-eps*gamma, 0., 0., -eps*gamma]
tc['b'] = lambda x,y : [eps]
tc['u'] = lambda x,y : [0.]
tc['f'] = lambda x,y : [0.5*(x**2+y**2)]
tc['bc_neumann'] = bc_neumann
# ...

# ...
PDE = basicPDE(geometry=geo, testcase=tc)
PDE.meanConstraint = False
# ...

# ...
V = PDE.space
V.nderiv_pts = 2
# ...

# ...
U   = PDE.unknown
rhs = PDE.rhs
# ...

# ...
PDE.assembly()
# ...

#-----------------------------------
def F(U,x,y):
    # ...
    P = V.get_points()
    x    = P[0,0,:]
    xdu  = P[0,1,:]
    xdv  = P[0,2,:]
    xduu = P[0,3,:]
    xduv = P[0,4,:]
    xdvv = P[0,5,:]

    y    = P[1,0,:]
    ydu  = P[1,1,:]
    ydv  = P[1,2,:]
    yduu = P[1,3,:]
    yduv = P[1,4,:]
    ydvv = P[1,5,:]

    jac = xdu * ydv - xdv * ydu
    # ...

    # ...
    D    = U.evaluate(patch_id=0, nderiv=2)

    _U   = D[0,0,:]
    Udu  = D[0,1,:]
    Udv  = D[0,2,:]
    Uduu = D[0,3,:]
    Uduv = D[0,4,:]
    Udvv = D[0,5,:]

    Udx =   ydv * Udu - ydu * Udv
    Udx /= jac
    Udy = - xdv * Udu + xdu * Udv
    Udy /= jac


    C1 = Uduu - xduu * Udx - yduu * Udy
    C2 = Uduv - xduv * Udx - yduv * Udy
    C3 = Udvv - xdvv * Udx - ydvv * Udy
    Udxx =   C1 * ydv**2    - 2 * C2 * ydu * ydv + C3 * ydu**2
    Udxx /= jac**2
    Udxy = - C1 * xdv * ydv + C2 *(xdu * ydv + xdv * ydu) - C3 * xdu * ydu
    Udxy /= jac**2
    Udyy =   C1 * xdv**2    - 2 * C2 * xdu * xdv + C3 * xdu**2
    Udyy /= jac**2
    # ...

    Hessian = Udxx * Udyy - Udxy**2
#    _F = sqrt ( abs(Hessian * rho1 (Udx,Udy) / (c_rho * rho0(x,y))) )
    _F = log ( abs(Hessian * rho1 (Udx,Udy) / (c_rho * rho0(x,y))) )

#    f_values = c_rho * rho0(x,y) / rho1 (Udx,Udy)
#    _F = - np.sqrt ( Udxx**2 + Udyy**2 + 2 * Udxy**2 + 2 * f_values )

    return [_F]
#-----------------------------------

#-----------------------------------
def plotMesh(PDE, ntx=60, nty=60):
    from matplotlib import pylab as plt

    geo = PDE.geometry
    patch_id = 0
    nrb   = geo[patch_id]

    C = np.zeros_like(nrb.points)

    _C = U.tomatrix(patch_id)
    shape = list(nrb.shape)
    C = np.zeros(shape+[3])
    C[...,0] = _C
    srf = cad_nurbs(nrb.knots, C, weights= nrb.weights)

    ub = srf.knots[0][0]
    ue = srf.knots[0][-1]
    vb = srf.knots[1][0]
    ve = srf.knots[1][-1]

    tx = np.linspace(ub, ue, ntx)
    ty = np.linspace(vb, ve, nty)

    nderiv = 1
    nderiv = 2

    # ...
    P    = nrb.evaluate_deriv(tx,ty,nderiv=nderiv)
    x    = P[0,:,:,0]
    xdu  = P[1,:,:,0]
    xdv  = P[2,:,:,0]
    xduu = P[3,:,:,0]
    xduv = P[4,:,:,0]
    xdvv = P[5,:,:,0]

    y    = P[0,:,:,1]
    ydu  = P[1,:,:,1]
    ydv  = P[2,:,:,1]
    yduu = P[3,:,:,1]
    yduv = P[4,:,:,1]
    ydvv = P[5,:,:,1]

    jac = xdu * ydv - xdv * ydu
    # ...

    # ...
    D    = srf.evaluate_deriv(tx,ty,nderiv=nderiv)
    Udu  = D[1,...,0]
    Udv  = D[2,...,0]
    Uduu = D[3,...,0]
    Uduv = D[4,...,0]
    Udvv = D[5,...,0]

    Udx =   ydv * Udu - ydu * Udv
    Udx /= jac
    Udy = - xdv * Udu + xdu * Udv
    Udy /= jac

    C1 = Uduu - xduu * Udx - yduu * Udy
    C2 = Uduv - xduv * Udx - yduv * Udy
    C3 = Udvv - xdvv * Udx - ydvv * Udy
    Udxx =   C1 * ydv**2    - 2 * C2 * ydu * ydv + C3 * ydu**2
    Udxy = - C1 * xdv * ydv + C2 *(xdu * ydv + xdv * ydu) - C3 * xdu * ydu
    Udyy =   C1 * xdv**2    - 2 * C2 * xdu * xdv + C3 * xdu**2
    # ...


    # ...
    fig = plt.figure()

#        Udx[:,0] = 0.
#        Udx[:,-1] = 1.
#        Udy[0,:] = 0.
#        Udy[-1,:] = 1.

    for i,v in enumerate(ty):
    #    phidx = Udu[:,i]
    #    phidy = Udv[:,i]

        phidx = Udx[:,i]
        phidy = Udy[:,i]

        plt.plot(phidx, phidy, '-b')

    for i,u in enumerate(tx):
    #    phidx = Udu[i,:]
    #    phidy = Udv[i,:]

        phidx = Udx[i,:]
        phidy = Udy[i,:]

        plt.plot(phidx, phidy, '-b')

    plt.show()
#-----------------------------------

# ...
def rhs_func(x,y):
    return F(U,x,y)
rhs.set_func(rhs_func)
# ...

# ...
u0 = lambda x,y : 0.5*(x**2+y**2)
PDE.interpolate(u0)

#PDE.plot(); plt.colorbar() ; plt.show()
# ...

# ...
list_Err = [1.e6]
t = 0.
i = 0
while (list_Err[-1] > rtol) and (i < maxiter):
    t += dt
    un = U.get()

    PDE.update()
    PDE.solve(rhs)
    dn = U.get()
    uh = un + dt * dn
    U.set(uh)

    err = np.linalg.norm(dn) / np.linalg.norm(un)
    list_Err.append(err)

    if verbose:
        print(i, ": ","   |F(x)| = ", list_Err[-1])
    i += 1
# ...

# ...
list_Err = np.asarray(list_Err[1:])
# ...

# ...
plotMesh(PDE, ntx=60, nty=60)
# ...

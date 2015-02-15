# -*- coding: UTF-8 -*-
#!/usr/bin/env python

import numpy as np
from sympy import *
from sympy.matrices import *
from scipy.optimize import root
from matplotlib import pyplot as plt
from scipy.integrate import quad

#from poisson_nonlin import poisson_picard as PDE_picard

#__all__ = ['genPoints', 'genFigure', 'genDomain', 'picard', 'picardTwoGrids', 'testcase']
#__all__ = ['genPoints', 'genFigure', 'genDomain', 'picard', 'picardTwoGrids', 'testcase']

sqrt = np.sqrt
abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)

# ------------------------------------------------------
# ... generate the boundary in the clock sens
def genPoints(R0=1., eps=0.32, k=1.7, d=0.33, mbnd=1500, m=500 \
                # needed for the newton algo to converge
                , rmin=0.67, rmax=1.33, delta=0.05 \
                , dd = 1.e-3 \
                , PLOT=False \
               ):
    print ("eps, k, d = ", eps, k, d)

    if (R0==1.) and (eps==0.32):
        d1, d2, d3 = [0.0753850296600659, -0.206294962187880, -0.0314337072805334]
    else:
        d1, d2, d3 = compute_ds(eps, k, d)

    print ("d1, d2, d3 = ", d1, d2, d3)

    psi   = lambda r,z: r**4/8 + d1 + d2 * r**2 + d3 * (r**4 - 4 * (r*z)**2)
    psidr = lambda r,z: 2*d2*r + d3*(4*r**3 - 8*r*z**2) + r**3/2
    psidz = lambda r,z: -8*d3*r**2*z

    # .....................................
    rgrid  = list(np.linspace(rmin, rmin+delta, mbnd)[:-1])
    rgrid += list(np.linspace(rmin+delta, rmax-delta, m))
    rgrid += list(np.linspace(rmax-delta, rmax, mbnd)[1:])
    rgrid = np.array(rgrid)

    zgrid = np.zeros_like(rgrid)

    # ...
    from pigasus.utils.impeqpy import impeqpy
    import pigasus.utils.impeqpy as impe
#    print "============================"
#    print impe.__file__
#    print "============================"
    level = 0.0
    imp=impeqpy(tol=1.e-9, maxniter = 300, verbose=False)
    imp.solve2Dx(psi,psidz,level,rgrid,zgrid)
    list_r = [R0-eps] ; list_z = [0.]
    for (r,z) in zip(rgrid, zgrid):
        if (not np.isnan(r)) and (not np.isnan(z)):
            list_r.append(r) ; list_z.append(z)
    list_r.append(R0+eps) ; list_z.append(0.)
    # ...
    # .....................................

    # .....................................
#    def Z2(R):
#        v  = 0.
#        v += d1/(4.*d3) * 1./R**2
#        v += d2/(4.*d3)
#        v += (1./8 + d3) / (4.*d3) * R**2
#        return v
#
#    def Z_plus(R):
#        return np.sqrt(Z2(R))
#
#    def Z_minus(R):
#        return - np.sqrt(Z2(R))
#
#    sol = root(Z2, rmin, jac=False)
#    Rmin = sol.x
#
#    sol = root(Z2, rmax, jac=False)
#    Rmax = sol.x
#
#    def dZdR(R):
#        gamma = (1./8 + d3) / (4.*d3)
#        alpha = d1/(4.*d3)
#        Z = Z_plus(R)
#        v = gamma * R - alpha / R**3
#        v /= Z
#        return v
#
#    def measure(R):
#        meas  = dZdR(R)**2
#    #    meas += 1.
#        return meas
#
#    def density(ti,tj):
#        return quad(measure, ti, tj)
#
#    def adaptive_mesh(n, xmin, xmax, amin, amax):
#        R = np.linspace(xmin, xmax, n)
#        D = []
#        for a,b in zip(R[:-1], R[1:]):
#            D.append(density(a,b)[0])
#
#        D = np.asarray(D)
#        m_total = D.sum()
#
#        M = np.zeros(n-1)
#        for i in range(0,n-1):
#            v = D[0:i]
#            M[i] = v.sum() / m_total
#
#        Rnew = (amax-amin)*M+amin
#        return Rnew
#
#    R = []
#    R += list(adaptive_mesh(mbnd, Rmin+dd, Rmin*(1.+delta), Rmin, Rmin*(1.+delta)))
#    R += list(adaptive_mesh(m, Rmin*(1.+delta), Rmax*(1.-delta), Rmin*(1.+delta), Rmax*(1.-delta)))
#    R += list(adaptive_mesh(mbnd, Rmax*(1.-delta), Rmax-dd, Rmax*(1.-delta), Rmax))
#    R = np.array(R)
#
#    Z = Z_plus(R)
#    R = np.array([Rmin] + list(R) + [Rmax])
#    Z = np.array([  0.] + list(Z) + [  0.])
#
#    list_r = R
#    list_z = Z
    # .....................................

    # ... y < 0 part
    rgrid = np.array(list_r); zgrid = np.array(list_z)

    _rgrid =  rgrid[::-1]
    _zgrid = -zgrid[::-1]
    # ...

    # ...
    if PLOT:
        import matplotlib.pyplot as plt
        plt.plot(rgrid, zgrid, '.b')
        plt.plot(_rgrid, _zgrid, '.k')

        tx = np.linspace(0.6,1.4, 300)
        ty = np.linspace(-0.6,0.6, 300)
        x,y = np.meshgrid(tx,ty)
        u = psi(x,y)
        levels = np.linspace(-0.04, 0.0, 100)
        CS = plt.contourf(x,y,u, levels)
        plt.colorbar()

        plt.show()

        genFigure(d=[d1,d2,d3], origin="upper")
    # ...

    r = list(rgrid) + list(_rgrid)
    z = list(zgrid) + list(_zgrid)

    return r,z
# ...
# ------------------------------------------------------


# ------------------------------------------------------
def genFigure(d=None, origin="upper"):
    #origin = 'lower'

    if d is None:
        d1, d2, d3 = [ 0.07538503, -0.20629496, -0.03143371]
    else:
        d1,d2,d3=d[0:3]

    import matplotlib.pyplot as plt

    # ITER and ASDEX-Upgrade
    tx = np.linspace(0.6,1.4, 300)
    ty = np.linspace(-0.6,0.6, 300)
    levels = np.linspace(-0.04, 0.0, 100)

#    # JET
#    levels = np.linspace(-0.045, 0., 100)
#    tx = np.linspace(0.6,1.4, 300)
#    ty = np.linspace(-0.65,0.65, 300)


    x,y = np.meshgrid(tx,ty)

    psi = lambda r,z: r**4/8 + d1 + d2 * r**2 + d3 * (r**4 - 4 * (r*z)**2)
    u = psi(x,y)
    #plt.contourf(x,y,u) ; plt.colorbar(); plt.show()



    CS = plt.contourf(x,y,u, levels
    #                        , colors = ('r', 'g', 'b') \
    #                        , origin=origin \
    #                        , extend='both' \
                      )
    plt.colorbar()

    CS2 = plt.contour(CS, levels=CS.levels[::10] \
                     , colors = 'k' \
                     , origin=origin \
                     , hold='on' \
                     , linewidths = (1,) \
                     )

    plt.show()
#    plt.pcolor(x,y,u, vmin=-0.04, vmax=0.01) ; plt.colorbar() ; plt.show()
# ------------------------------------------------------

# ------------------------------------------------------
class genDomain(object):
    def __init__(self, R0=1, eps=0.32, mbnd=1500, m=500 \
                        # needed for the newton algo to converge
                        , rmin=0.67, rmax=1.33, delta=0.05 \
                        , PLOT=False):

        r,z = genPoints(R0=R0, eps=eps, mbnd=mbnd, m=m \
                        , rmin=rmin, rmax=rmax, delta=delta \
                        , PLOT=PLOT\
                       )
        self.boundary = [r,z]
# ------------------------------------------------------
def compute_ds(epsilon, kappa, delta):
    # ... OLD VERSION
#    from scipy import matrix
#    from scipy.linalg import inv
#    M = np.zeros((3,3))
#    M[:,0] = 1.
#    M[0,1] = (1+eps)**2
#    M[1,1] = (1-eps)**2
#    M[2,1] = (1-d*eps)**2
#    M[0,2] = (1+eps)**4
#    M[1,2] = (1-eps)**4
#    M[2,2] = (1-d*eps)**4 - 4 * ( (1-d*eps)*k*eps )**2
#    Y = np.zeros(3)
#    Y[0] = -(1./8) * (1+eps)**4
#    Y[1] = -(1./8) * (1-eps)**4
#    Y[2] = -(1./8) * (1-d*eps)**4
#    A = matrix(M)
#    invA = inv(A)
#    X = invA.dot(Y)
    # ...
    def compute_M():
        e = Symbol('e')
        k = Symbol('k')
        d = Symbol('d')

        A = Matrix([  [1,  (1+e)**2,                         (1+e)**4] \
                    , [1,  (1-e)**2,                         (1-e)**4] \
                    , [1,(1-d*e)**2,(1-d*e)**4 - 4.*((1.-d*e)*k*e)**2] \
                   ])

        Ainv =  A.inv()
        M = lambdify((e,k,d), Ainv)
        return M

    M = compute_M()

    Y = np.zeros(3)
    Y[0] = -(1./8) * (1+epsilon)**4
    Y[1] = -(1./8) * (1-epsilon)**4
    Y[2] = -(1./8) * (1-delta * epsilon)**4

    D= M(epsilon, kappa, delta).dot(Y)
    d1 = D[0,0]
    d2 = D[0,1]
    d3 = D[0,2]
    return d1, d2, d3

class testcase(object):
    def __init__(self, TEST):
        initTEST = getattr(self, 'initTEST%d' % TEST)
        initTEST()

    def initTEST1(self):
        """
        ITER relevant parameters
        d1, d2, d3 = [ 0.0753850296600659 -0.206294962187880 -0.0314337072805334]
        """
        d1, d2, d3 = compute_ds(eps=0.32, k=1.7, d=0.33)

        # ...
        F = lambda psi,x,y : x**2
        psi = lambda r,z: r**4/8 + d1 + d2 * r**2 + d3 * (r**4 - 4 * (r*z)**2)
        # ...

        self.F = F

#if __name__ == '__main__':
#    from caid.cad_geometry import square
#    from matplotlib import pylab as plt


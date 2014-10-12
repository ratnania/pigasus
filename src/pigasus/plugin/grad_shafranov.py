# -*- coding: UTF-8 -*-
#!/usr/bin/env python

import numpy as np
#from poisson_nonlin import poisson_picard as PDE_picard

#__all__ = ['genPoints', 'genFigure', 'genDomain', 'picard', 'picardTwoGrids', 'testcase']
#__all__ = ['genPoints', 'genFigure', 'genDomain', 'picard', 'picardTwoGrids', 'testcase']

sqrt = np.sqrt
abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)

# ------------------------------------------------------
# ... generate the boundary in the clock sens
def genPoints(R0=1, eps=0.32, k=1.7, d=0.33, mbnd=1500, m=500 \
                # needed for the newton algo to converge
                , rmin=0.67, rmax=1.33, delta=0.05 \
                , PLOT=False \
               ):

    if (R0==1.) and (eps==0.32):
        d1, d2, d3 = [ 0.07538503, -0.20629496, -0.03143371]
    else:
        d1, d2, d3 = compute_ds(eps=eps, k=k, d=d)

    psi   = lambda r,z: r**4/8 + d1 + d2 * r**2 + d3 * (r**4 - 4 * (r*z)**2)
    psidr = lambda r,z: 2*d2*r + d3*(4*r**3 - 8*r*z**2) + r**3/2
    psidz = lambda r,z: -8*d3*r**2*z

    rgrid  = list(np.linspace(rmin, rmin+delta, mbnd)[:-1])
    rgrid += list(np.linspace(rmin+delta, rmax-delta, m))
    rgrid += list(np.linspace(rmax-delta, rmax, mbnd)[1:])
    rgrid = np.array(rgrid)

    zgrid = np.zeros_like(rgrid)

    # ...
    from pigasus.utils.impeqpy import impeqpy
    level = 0.0
    imp=impeqpy()
    imp.solve2Dx(psi,psidz,level,rgrid,zgrid)
    list_r = [R0-eps] ; list_z = [0.]
    for (r,z) in zip(rgrid, zgrid):
        if (not np.isnan(r)) and (not np.isnan(z)):
            list_r.append(r) ; list_z.append(z)
    list_r.append(R0+eps) ; list_z.append(0.)
    # ...

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

    tx = np.linspace(0.6,1.4, 300)
    ty = np.linspace(-0.6,0.6, 300)
    x,y = np.meshgrid(tx,ty)

    psi = lambda r,z: r**4/8 + d1 + d2 * r**2 + d3 * (r**4 - 4 * (r*z)**2)
    u = psi(x,y)
    #plt.contourf(x,y,u) ; plt.colorbar(); plt.show()

    #levels = np.linspace(-0.04, 0.01, 50)
    levels = np.linspace(-0.04, 0.0, 100)
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
def compute_ds(eps=0.32, k=1.7, d=0.33):
    from scipy import matrix
    from scipy.linalg import inv
    M = np.zeros((3,3))
    M[:,0] = 1.
    M[0,1] = (1+eps)**2
    M[1,1] = (1-eps)**2
    M[2,1] = (1-d*eps)**2
    M[0,2] = (1+eps)**4
    M[1,2] = (1-eps)**4
    M[2,2] = (1-d*eps)**4 - 4 * ( (1-d*eps)*k*eps )**2
    Y = np.zeros(3)
    Y[0] = -(1./8) * (1+eps)**4
    Y[1] = -(1./8) * (1-eps)**4
    Y[2] = -(1./8) * (1-d*eps)**4
    A = matrix(M)
    invA = inv(A)
    X = invA.dot(Y)
    return X


class testcase(object):
    def __init__(self, TEST):
        initTEST = getattr(self, 'initTEST%d' % TEST)
        initTEST()

    def initTEST1(self):
        """
        ITER relevant parameters
        d1, d2, d3 = [ 0.07538503, -0.20629496, -0.03143371]
        """
        d1, d2, d3 = compute_ds(eps=0.32, k=1.7, d=0.33)

        # ...
        F = lambda psi,x,y : x**2
        psi = lambda r,z: r**4/8 + d1 + d2 * r**2 + d3 * (r**4 - 4 * (r*z)**2)
        # ...

        self.F = F

#if __name__ == '__main__':
#    from igakit.cad_geometry import square
#    from matplotlib import pylab as plt


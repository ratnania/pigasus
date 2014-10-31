# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from scipy.sparse.linalg import spsolve
from pigasus.gallery.poisson_nonlin import poisson_picard as PDE_picard
#from basicPDE_nonlin import basicPDE_picard as PDE_picard
from pigasus.gallery.poisson import *
from time import time
from caid.cad_geometry import cad_geometry, cad_nurbs

__all__ = ['picard', 'picardTwoGrids', 'testcase']

sqrt = np.sqrt
abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)

class picard(PDE_picard):
    """
    A multidimentional nonlinear Monge Ampere class solver using Picard algorithm.
        >>> import caid.cad_geometry  as cg
        >>> from caid.cad_geometry import line
        >>> import pylab                as pl


    """

    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, *args, **kwargs):
        """Creates a nonlinear poisson PDE solver based on Picard algorithm.

        geometry:
            The geometry must be an object cad_geometry.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.gallery.poisson`.

        """

        # ...
        PDE_picard.__init__(self, *args, **kwargs)
        # ...

        # ...
        V = self.space
        V.nderiv_pts = 2
        # ...
    #-----------------------------------

    #-----------------------------------
    def initialize(self, u0=None):
        U = self.unknown
        if u0 is None:
            U.set(np.zeros(U.size))
            return
        if u0.__class__.__name__=="ndarray":
            U.set(u0)
        else:
            self.interpolate(u0, field=U)
    #-----------------------------------

    #-----------------------------------
    def solve(self, rho0, rho1, c_rho=None, u0=None, maxiter=100, rtol=1.e-6, rtol2=1.e-6 \
              , verbose=False, update=False):
        """
        solves the nonlinear poisson equation using PIcard algorithm
        rho0:
            the initial density
        rho1:
            the new density
        u0:
            this is the initial value for u. Default: all B-splines coeff = 0
        maxiter:
            the maximum number of iterations for the Picard algorithm. Default 100
        rtol:
            the relative tolerance. Default 1.e-6
        verbose:
            True => print the error for each iteration

        Returns:
            The residual error (as a numpy array)

        """
        PDE = self
        V   = self.space

#        print ">> solve-Monge Ampere"

        # ... compute the ratio int rho1 / int rho0
        if c_rho is None:
            # assembly the stifness matrix and bc terms
#            print "Entering assembly"
            PDE_picard.assembly(self, update=update)
#            print "Leaving assembly"
            if self.Dirichlet:
                U = self.unknown_dirichlet
            else:
                U = self.unknown
            U.reset()
            u_sqr = lambda x,y : [sqrt(rho0(x,y))]
            C0 = 1./PDE.norm(exact=u_sqr)**2
            u_sqr = lambda x,y : [sqrt(rho1(x,y))]
            C1 = 1./PDE.norm(exact=u_sqr)**2
            c_rho = C0/C1
        # ...

        self.rho0  = rho0
        self.rho1  = rho1
        self.c_rho = c_rho

        # ...
        from pigasus.fem.utils import function
        def _F (U,x,y):
            D    = U.evaluate(nderiv=2, parametric=False)

            _U   = D[0,0,:]
            Udx  = D[0,1,:]
            Udy  = D[0,2,:]
            Udxx = D[0,3,:]
            Udxy = D[0,4,:]
            Udyy = D[0,5,:]

            f_values = c_rho * rho0(x,y) / rho1 (Udx,Udy)
            return [-sqrt ( Udxx**2 + Udyy**2 + 2 * Udxy**2 + 2 * f_values)]
        # ...
        func_F = function(_F, fields=[U])

#        # ...
#        def F(U,x,y):
#
#            # ...
#            D    = U.evaluate(nderiv=2, parametric=False)
#
#            _U   = D[0,0,:]
#            Udx  = D[0,1,:]
#            Udy  = D[0,2,:]
#            Udxx = D[0,3,:]
#            Udxy = D[0,4,:]
#            Udyy = D[0,5,:]
#
#            f_values = c_rho * rho0(x,y) / rho1 (Udx,Udy)
#            _F = - np.sqrt ( Udxx**2 + Udyy**2 + 2 * Udxy**2 + 2 * f_values )
#
#            return [_F]
#        # ...

        return PDE_picard.solve(  self, func_F \
                             , u0=u0, maxiter=maxiter, rtol=rtol, rtol2=rtol2 \
                             , verbose=verbose, update=update)
    #-----------------------------------

    #-----------------------------------
    def plotMesh(self, ntx=60, nty=60):
        from matplotlib import pylab as plt

        geo = self.geometry
        patch_id = 0
        nrb   = geo[patch_id]

        C = np.zeros_like(nrb.points)
        if self.Dirichlet:
            U = self.unknown_dirichlet
        else:
            U = self.unknown
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

#        plt.show()
    #-----------------------------------

    #-----------------------------------
    def transferSolution(self, geo_H, U_H, geo_h, U_h):
        patch_id = 0
        nrb_H = geo_H[patch_id]
        nrb_h = geo_h[patch_id]

        # degrees of the input surface
        px_in = nrb_H.degree[0]
        py_in = nrb_H.degree[1]

        # degrees of the output surface
        px = nrb_h.degree[0]
        py = nrb_h.degree[1]

        C = np.zeros_like(nrb_H.points)
        shape = list(nrb_H.shape)
        C = np.zeros(shape+[3])
        _C = U_H.tomatrix(patch_id)
        C[...,0] = _C
        srf = cad_nurbs(nrb_H.knots, C, weights= nrb_H.weights)

        geo_f = cad_geometry()
        geo_f.append(srf)

        list_t = []
        for axis in range(0, nrb_H.dim):
            tx = [t for t in nrb_h.knots[axis] if t not in nrb_H.knots[axis]]
            list_t.append(tx)

        geo_f.refine(list_t = list_t, list_p=[px-px_in, py-py_in])
        uh = geo_f[0].points[...,0]
        U_h.frommatrix(patch_id, uh)
    #-----------------------------------

    #-----------------------------------
    def Err_func(self, x,y):
        V = self.space
        if self.Dirichlet:
            U = self.unknown_dirichlet
        else:
            U = self.unknown
        rho1 = self.rho1
        c_rho = self.c_rho

        # ...
        D    = U.evaluate(nderiv=2, parametric=False)
        U.reset()

        _U   = D[0,0,:]
        Udx  = D[0,1,:]
        Udy  = D[0,2,:]
        Udxx = D[0,3,:]
        Udxy = D[0,4,:]
        Udyy = D[0,5,:]
        # ...

    #    _F = c_rho * rho0(x,y) / rho1 (udx,udy) - ( udxx*udyy - udxy**2 )
        _F = ( Udxx*Udyy - Udxy**2 ) * rho1 (Udx,Udy) / c_rho

        return [_F]
    #-----------------------------------

class picardTwoGrids(object):
    """
    A multidimentional nonlinear Monge Ampere class solver using a two grids Picard algorithm.
        >>> import caid.cad_geometry  as cg
        >>> from caid.cad_geometry import line
        >>> import pylab                as pl


    """
    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, geometry, geometry_H, bc_dirichlet=None, bc_neumann=None,
                 AllDirichlet=None, Dirichlet=None, metric=None, solverInfo=None):
        """Creates a nonlinear poisson PDE solver based on Picard algorithm.

        geometry:
            The geometry must be an object cad_geometry.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.gallery.poisson`.

        """

        # ...
        self.PDE_h = picard( geometry \
                        , bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann \
                        , AllDirichlet=AllDirichlet, Dirichlet=Dirichlet \
                        , metric=metric, solverInfo=solverInfo)


        self.PDE_H = picard( geometry_H \
                        , bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann \
                        , AllDirichlet=AllDirichlet, Dirichlet=Dirichlet \
                        , metric=metric, solverInfo=solverInfo)
        # ...

    #-----------------------------------
    def initialize(self, u0=None):
        U = self.PDE_H.unknown
        if u0 is None:
            U.set(np.zeros(U.size))
            return
        if u0.__class__.__name__=="ndarray":
            U.set(u0)
        else:
            self.PDE_H.interpolate(u0, field=U)
    #-----------------------------------

    #-----------------------------------
    def solve(self, rho0, rho1, c_rho=None, u0=None \
              , maxiter=[50,100], rtol=[1.e-5, 1.e-6], rtol2=[1.e-5,1.e-6] \
              , verbose=0, update=False):

        PDE_h = self.PDE_h
        PDE_H = self.PDE_H

        geo_h = PDE_h.geometry
        geo_H = PDE_H.geometry

        # ...
        Dirichlet = PDE_h.Dirichlet
        if Dirichlet:
            U_h = PDE_h.unknown_dirichlet
            U_H = PDE_H.unknown_dirichlet
        else:
            U_h = PDE_h.unknown
            U_H = PDE_H.unknown
        # ...

        # ...
        if verbose > 0:
            print("*****************************")
            tb = time()

        Errors_H, ErrorsH1_H = PDE_H.solve(  rho0, rho1 \
                    , c_rho=c_rho, u0=u0 \
                    , maxiter=maxiter[0], rtol=rtol[0], rtol2=rtol2[0] \
                    , verbose=(verbose==2), update=update)

        if verbose > 0:
            te = time()
            print("Coarse solver converges after ", len(Errors_H) \
                    , " with final error ", Errors_H[-1] \
                    , " with final H1-error ", ErrorsH1_H[-1])
            print("Elapsed time ", te-tb)
            print("*****************************")
        # ...

        # ...
        PDE_H.transferSolution(geo_H, U_H, geo_h, U_h)
        u0 = U_h.get()
        # ...

        # ...
        if verbose > 0:
            print("*****************************")
            tb = time()

        Errors_h, ErrorsH1_h = PDE_h.solve(  rho0, rho1 \
                                           , c_rho=c_rho, u0=u0 \
                                           , maxiter=maxiter[1] \
                                           , rtol=rtol[1], rtol2=rtol2[1] \
                                           , verbose=(verbose==1))

        if verbose > 0:
            te = time()
            print("Monge-Ampere eq. converges after ", len(Errors_h) \
                    , " with final error ", Errors_h[-1] \
                    , " with final H1-error ", ErrorsH1_h[-1])
            print("Elapsed time ", te-tb)
            print("*****************************")
        # ...

        return Errors_h, ErrorsH1_h, Errors_H, ErrorsH1_H

    #-----------------------------------
    def plotMesh(self, ntx=60, nty=60):
        self.PDE_h.plotMesh(ntx=ntx, nty=nty)
    #-----------------------------------

    #-----------------------------------
    def Err_func(self, x,y):
        self.PDE_h.Err_func(x,y)
    #-----------------------------------

    #-----------------------------------
    def free(self):
        self.PDE_h.free()
        self.PDE_H.free()
    #-----------------------------------

class testcase(object):
    def __init__(self, TEST):
        initTEST = getattr(self, 'initTEST%d' % TEST)
        initTEST()

    def initTEST1(self):
        C0 = 1.0

        # ... test 1
        rho0 = lambda x,y : 1.
        C1   = 0.616805883732
        t = 0.5
        rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-t)**2+(y-0.5)**2-0.09)))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST2(self):
        C0 = 1.0

        # ... test 2
        rho0 = lambda x,y : 1.
        C1   =  1.75484181939
        rho1 = lambda x,y : ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST3(self):
        C0 = 1.0

        # ... test 3
        rho0 = lambda x,y : 1.
        C1   = 0.285547502263
        rho1 = lambda x,y : (1. + 10*exp(-50*(y-0.5-0.25*sin(2*pi*x))**2))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST4(self):
        C0 = 1.0

        # ...
        rho0 = lambda x,y : 1.
        t = 0.25
        C1 = 0.702563292151
        rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-0.25*cos(2*pi*t))**2 \
                            - (y-0.5-0.5 *sin(2*pi*t))**2 \
                            - 0.01) ))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST5(self):
        C0 = 1.0

        # ...
        rho0 = lambda x,y : 1.
        t = 1./3
        C1 = 0.831806957866
        rho1 = lambda x,y : ( 1. + 5*exp(-50*abs(y-0.5-0.25*sin(2*pi*x)*sin(2*pi*t))))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho


    def initTEST6(self):
        C0 = 1.0

        # ...
        rho0 = lambda x,y : 1.
        gamma = 5.
        lamb = 100.
        t = 0.75
        C1 = 0.832943327557
        x0 = t ; y0 = 0.2 + 0.5 * t ; x1 = 1. - t ; y1 = 0.8 - 0.5 * t
        u0 = lambda x,y : gamma * sech(lamb * ( x - x0 + y - y0 ))
        u1 = lambda x,y : gamma * sech(lamb * ( x - x1 + y - y1 ))
        rho1 = lambda x,y : ( 1. + u0(x,y) + u1(x,y))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST7(self):
        C0 = 1.0

        # ... test7
        xc = 0.7 ; yc = 0.5
        C1 = 0.281648379406

        rho0 = lambda x,y : 1.
        r = lambda s,t : sqrt( (s-xc)**2 + (t-yc)**2 )
        theta = lambda s,t : atan(t-yc,s-xc)
        def rho1(s,t):
            r_ = r(s,t) ;  t_ = theta(s,t)
            val = C1 * (1. + 9./(1. + (10*r_*cos(t_-20*r_**2))**2) )
            return val
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST8(self):
        C0 = 1.0

        # ... circle
        rho0 = lambda x,y : 1./pi
        C1   = 0.227475185932
        rho1 = lambda x,y : C1 * (1. + 5*exp(-25*abs((x-0.)**2+(y-0.)**2-0.4)))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST9(self):
        C0 = 1.0

        # ... quart_circle
        rho0 = lambda x,y : 4./pi
        C1   = 2.91639889933
        rho1 = lambda x,y : C1 * ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

    def initTEST10(self):
        C0 = 1.0

        # ... annulus
        C0   = 0.424413181542
        rho0 = lambda x,y : C0 * 1.
        C1   = 0.733393862165
        rho1 = lambda x,y : C1 * ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
        # ...

        c_rho = C0/C1

        self.rho0 = rho0
        self.rho1 = rho1
        self.c_rho = c_rho

#    def initTEST(self):
#        C0 = 1.0
#
#
#
#        c_rho = C0/C1
#
#        self.rho0 = rho0
#        self.rho1 = rho1
#        self.c_rho = c_rho


#if __name__ == '__main__':
#    from caid.cad_geometry import square
#    from matplotlib import pylab as plt
#    from time import time
#
#    abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
#    pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
#    sech = lambda x: 1./cosh(x)
#
#    verbose = False
#    withTwoGrids = True
#
##    p_H = [ 5, 5]
##    p_h = [ 5, 5]
#
#    p_H = [ 3, 3]
#    p_h = [ 3, 3]
#
##    p_H = [ 2, 2]
##    p_h = [ 2, 2]
#
#    # TEST 1
##    # p = 5
##    rtol_H = 1.e-4
##    rtol_h = 1.e-8
#
#    # p = 3
##    rtol_H = 1.e-4
##    rtol_h = 1.e-6
#
#    # p = 2
##    rtol_H = 1.e-4
##    rtol_h = 1.e-4
#
#    # TEST 3
#    # p = 3, 5
#    rtol_H = 1.e-3
#    rtol_h = 1.e-3
#
##    # p = 2
###    rtol_H = 1.e-2
###    rtol_h = 1.e-2
#
#
#    rtol2_H = 1.e-6
##    rtol2_h = 1.e-6
#    rtol2_h = 1.e-9
#
#    maxiter_H = 40
#    maxiter_h = 40
#
#    n_H = [7,7]
#    nstage =  1
#    n_h = []
#    for axis in range(0,2):
#        n = n_H[axis]
#        for i in range(0, nstage):
#            n = 2*n+1
#        n_h.append(n)
#
#    if withTwoGrids:
#        print ">>>> coarse grid ", n_H, " with splines of degree ", p_H
#    print ">>>> fine   grid ", n_h, " with splines of degree ", p_h
#
#    if withTwoGrids:
#        geo_H = square(n=n_H, p=p_H)
#    geo_h = square(n=n_h, p=p_h)
#    # ...
#
#    # ...
#    C0   = 1.0
#    rho0 = lambda x,y : 1.
#
##    C1   = 0.616805883732
##    t = 0.5
##    rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-t)**2+(y-0.5)**2-0.09)))
#
##    C1   =  1.75484181939
##    rho1 = lambda x,y : ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
#
#    # ... test7
#    xc = 0.7 ; yc = 0.5
#    C1 = 0.281648379406
#
#    r = lambda s,t : sqrt( (s-xc)**2 + (t-yc)**2 )
#    theta = lambda s,t : atan(t-yc,s-xc)
#    def rho1(s,t):
#        r_ = r(s,t) ;  t_ = theta(s,t)
#        val = C1 * (1. + 9./(1. + (10*r_*cos(t_-20*r_**2))**2) )
#        return val
#    # ...
#
#    c_rho = C0/C1
#    # ...
#
#    # ...
#    # values of gradu.n at the boundary
#    # ...
#    def func_g(x,y):
#        return [x,y]
#    # ...
#
#    # ...
#    # values of u at the boundary
#    # ...
#    bc_neumann={}
#    bc_neumann [0,0] = func_g
#    bc_neumann [0,1] = func_g
#    bc_neumann [0,2] = func_g
#    bc_neumann [0,3] = func_g
#    # ...
#
#    PDE_h = picardTwoGrids(geometry=geo_h, geometry_H=geo_H, bc_neumann=bc_neumann)
#
#    Errors_h, ErrorsH1_h, Errors_H, ErrorsH1_H = \
#            PDE_h.solve(  rho0, rho1, c_rho=None, u0=None \
#                        , maxiter=[maxiter_H,maxiter_h] \
#                        , rtol=[rtol_H,rtol_h] \
#                        , rtol2=[rtol2_H,rtol2_h] \
#                        , verbose=1)
#    # ...
##    tc = {}
##    tc['A'] = lambda x,y : [1., 0., 0., 1.]
##    tc['b'] = lambda x,y : [1.e-3]
##    tc['u'] = lambda x,y : [0.]
##    tc['f'] = lambda x,y : [0.]
##    tc['bc_neumann'] = bc_neumann
#    # ...
#
##    PDE_H = picard(geometry=geo_H, testcase=tc)
##    PDE_h = picard(geometry=geo_h, testcase=tc)
#
##    if withTwoGrids:
##        PDE_H = picard(geometry=geo_H, bc_neumann=bc_neumann)
##    PDE_h = picard(geometry=geo_h, bc_neumann=bc_neumann)
##
##    # ...
##    print ">>> Solving using Picard <<<"
##    # ...
##    if withTwoGrids:
##        if PDE_H.Dirichlet:
##            U_H = PDE_H.unknown_dirichlet
##        else:
##            U_H = PDE_H.unknown
##
##    if PDE_h.Dirichlet:
##        U_h = PDE_h.unknown_dirichlet
##    else:
##        U_h = PDE_h.unknown
##    # ...
##
##    # ...
##    if withTwoGrids:
##        print "*****************************"
##        tb = time()
##        Errors_H, ErrorsH1_H = PDE_H.solve(  rho0, rho1, c_rho=None, u0=None \
##                    , maxiter=maxiter_H, rtol=rtol_H, rtol2=rtol2_h, verbose=verbose)
##        te = time()
##        print "Coarse solver converges after ", len(Errors_H) \
##                , " with final error ", Errors_H[-1] \
##                , " with final H1-error ", ErrorsH1_H[-1]
##        print "Elapsed time ", te-tb
##        print "*****************************"
##
##        PDE_H.transferSolution(geo_H, U_H, geo_h, U_h)
##        u0 = U_h.get()
##    else:
##        u0 = np.zeros(PDE_h.size)
##
##    print "*****************************"
##    tb = time()
##    Errors_h, ErrorsH1_h = PDE_h.solve(  rho0, rho1, c_rho=None, u0=u0 \
##                , maxiter=maxiter_h, rtol=rtol_h, rtol2=rtol2_h, verbose=verbose)
##    te = time()
##    print "Monge-Ampere eq. converges after ", len(Errors_h) \
##            , " with final error ", Errors_h[-1] \
##            , " with final H1-error ", ErrorsH1_h[-1]
##    print "Elapsed time ", te-tb
##    print "*****************************"
##
##    if withTwoGrids:
##        uH = U_H.get()
##    uh = U_h.get()
##
##    if withTwoGrids:
##        print "Error-coarse        ", np.abs(1.-PDE_H.norm(exact=PDE_H.Err_func))
##    print "Error-fine          ", np.abs(1.-PDE_h.norm(exact=PDE_h.Err_func))
##
##    if withTwoGrids:
##        U_H.set(uH)
##    U_h.set(uh)
##    # ...
##
##    # ...
###    PDE_H.plotMesh(ntx=60, nty=60)
###    PDE_h.plotMesh(ntx=60, nty=60)
##    # ...
##
##    np.savetxt("Errors.txt", np.asarray(Errors_h))

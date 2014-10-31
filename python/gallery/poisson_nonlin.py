# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from scipy.sparse.linalg import spsolve
from .poisson import *
from pigasus.fem.basicPDE import *
from numpy import abs

__all__ = ['poisson_picard', 'poisson_newton']

class poisson_picard(poisson):
    """
    A multidimentional nonlinear Poisson class solver using Picard algorithm.
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

        poisson.__init__(self, *args, **kwargs)

        # ...
    #-----------------------------------

    #-----------------------------------
    def initialize(self, u0=None):
        U = self.unknown
        if u0 is None:
            U.set(np.zeros(U.size))
            return
#        self.project(u0, field=U)
        self.interpolate(u0, field=U)
    #-----------------------------------

    #-----------------------------------
    def assembly(self, f=None, update=False):
        poisson.assembly(self, f=f, update=update)
    #-----------------------------------

    #-----------------------------------
    def solve(self, F, u0=None, maxiter=100, rtol=1.e-6, rtol2=1.e-6 \
              , verbose=False, update=False):
        """
        solves the nonlinear poisson equation using PIcard algorithm
        F:
            the rhs. it can be any function F(U, gradU, ..., x,y)
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
        # assembly the stifness matrix and bc terms
        poisson.assembly(self, update=update)

        # project u0 onto the discrete vectorial space
        self.initialize(u0=u0)

        # ...
        PDE  = self
        V = PDE.space
        un   = PDE.unknown
        rhs  = self.rhs
        # ...

        rhs.func = F

        # ...
        from time import time
        list_Err   = [1.e6]
        list_ErrH1 = [1.e6]
        un_values = un.get()
        normH1_old = np.dot(PDE.dot(un.get()), un.get())
        i = 0
        if verbose:
            tb = time()
        while (list_Err[-1] > rtol) and (list_ErrH1[-1] > rtol2) and (i < maxiter):
            U_old_values = un.get()
#            print "-------"
#            print "solve"
#            import matplotlib.pyplot as plt
##            Phi = PDE.G_W
#            Phi = PDE.unknown_dirichlet
##            Phi.plot(withpcolor=True) ; plt.colorbar() ; plt.show()
#            Phi.fast_plot() ; plt.colorbar() ; plt.show()
#            print "-------"

            # assembly the right hand side
            rhs.reset()
            self.update()
            # solve and update unew
            poisson.solve(self, rhs)

            U_values = un.get()
            err = np.linalg.norm(U_values-U_old_values)
            list_Err.append(err)

            normH1 = np.dot(PDE.dot(un.get()), un.get())
            list_ErrH1.append(np.abs(normH1-normH1_old))

            normH1_old = normH1

            i += 1
            if verbose:
                print(i, ": ","   |F(x)| = ", list_Err[-1],"   |DF(x)| = ", list_ErrH1[-1])
        if verbose:
            te = time()
            print(">> Elapsed time ", te-tb)

        list_Err   = np.asarray(list_Err[1:])
        list_ErrH1 = np.asarray(list_ErrH1[1:])
        return list_Err, list_ErrH1
    #-----------------------------------

class poisson_newton(poisson):
    """
    A multidimentional nonlinear Poisson class solver using Picard algorithm.
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
        try:
            geometry = kwargs['geometry']
        except:
            pass
        # ...
        dim = geometry.dim
        if dim == 1:
            func_one   = lambda x : [ 1. ]
            func_zero  = lambda x : [ 0. ]
            func_stiff  = lambda x : [ 1. ]
        if dim == 2:
            func_one   = lambda x,y : [ 1. ]
            func_zero  = lambda x,y : [ 0. ]
            func_stiff  = lambda x,y : [  1., 0. \
                                        , 0., 1. ]
        if dim == 3:
            func_one   = lambda x,y,z : [ 1. ]
            func_zero  = lambda x,y,z : [ 0. ]
            func_stiff  = lambda x,y,z : [ 1., 0., 0. \
                                         , 0., 1., 0. \
                                         , 0., 0., 1. ]
        # ...

        # ...
        tc_d = {}
        tc_d['A'] = func_stiff
        tc_d['b'] = func_zero
        try:
            tc_d['AllDirichlet'] = kwargs['AllDirichlet']
        except:
            pass
        try:
            tc_d['bc_dirichlet'] = kwargs['bc_dirichlet']
        except:
            pass
        try:
            tc_d['bc_neumann'] = kwargs['bc_neumann']
        except:
            pass
        try:
            tc_d['Metric'] = kwargs['Metric']
        except:
            pass
        # ...

        # ...
        poisson.__init__(self, *args, **kwargs)
        self.Dn = basicPDE(geometry=geometry, testcase=tc_d)
        # ...

        # ...
    #-----------------------------------

#    #-----------------------------------
#    def __del__(self):
#        self.Dn.__del__()
#        poisson.__del__(self)
#    #-----------------------------------

    #-----------------------------------
    def free(self):
        self.Dn.free()
        poisson.free(self)
    #-----------------------------------

    #-----------------------------------
    def initialize(self, u0=None):
        U = self.unknown
        if u0 is None:
            U.set(np.zeros(U.size))
            return
#        self.project(u0, field=U)
        self.interpolate(u0, field=U)
    #-----------------------------------

    #-----------------------------------
    def solve(self, F, dF, u0=None, maxiter=100, rtol=1.e-6 \
              , verbose=False, update=False):
        """
        solves the nonlinear poisson equation using PIcard algorithm
        F:
            the rhs. it can be any function F(U, gradU, ..., x,y)
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
        # assembly the stifness matrix and bc terms
        poisson.assembly(self, update=update)
        self.Dn.assembly()

        # project u0 onto the discrete vectorial space
        self.initialize(u0=u0)

        En = self
        Dn = self.Dn

        # ...
        if En.Dirichlet:
            U = En.unknown_dirichlet
        else:
            U = En.unknown
        # ...

        # ...
        # current values
        un = En.unknown
        # unew-un
        dn = Dn.unknown
        # get the right hand side
        rhs = En.rhs
        # redefine the right hand side function
        def rhs_func(x,y):
            return F(U,x,y)
        rhs.set_func(rhs_func)
        def Mn_func(x,y):
            return dF(U,x,y)
        # get the mass operator
        Mn = Dn.mass
        # redefine the mass function
        Mn.set_func (Mn_func)
        # ...

        # ...
        from time import time
        dn.reset()
        list_Err    = [1.e6]
        list_ErrH1  = [1.e6]
        un_values = un.get()
        i = 0
        tb = time()
        while (list_Err[-1] > rtol) and (i < maxiter):
            # assembly the right hand side
            rhs.reset()
            En.update()
            Dn.assembly()
            # compute the right hand side
            g = rhs - En.dot (un)
            # solve and update unew
            Dn.solve (g)
            un += dn

            err = np.linalg.norm(dn.get())
            list_Err.append(err)

            err = np.dot(self.Dn.dot(dn.get()), dn.get())
            list_ErrH1.append(abs(err))

            i += 1
            if verbose:
                print(i, ": ","   |F(x)| = ", list_Err[-1],"   |DF(x)| = ", list_ErrH1[-1])
        te = time()
        print(">> Elapsed time ", te-tb)

        list_Err   = np.asarray(list_Err[1:])
        list_ErrH1 = np.asarray(list_ErrH1[1:])
        return list_Err, list_ErrH1
    #-----------------------------------

if __name__ == '__main__':
    from caid.cad_geometry import circle
    from matplotlib import pylab as plt

    sin = np.sin ; cos = np.cos ; exp = np.exp ; log = np.log ; sqrt = np.sqrt ; pi = np.pi

    nx = 15 ; ny = 15
    px = 2  ; py = 2

    geo = circle(radius=1./sqrt(2), n=[nx,ny], p=[px,py])

    # ...
    u_exact = lambda x,y : [- 2.0 * log ( x**2 + y**2 + 0.5 )]

    def F(U,x,y):
        _U = U.evaluate()
        return [4. * exp(_U)]

    def dF (U,x, y):
        _U = U.evaluate()
        return[-4 * exp(_U)]
    # ...

    AllDirichlet = True
    PDE_picard = poisson_picard(  geometry=geo \
                         , AllDirichlet=AllDirichlet )
    PDE_newton = poisson_newton(  geometry=geo \
                         , AllDirichlet=AllDirichlet )


    # ...
    print(">>> Solving using Picard <<<")
    # ...
    PDE = PDE_picard
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown
    # ...
    PDE_picard.solve(F, u0=None, maxiter=100, rtol=1.e-6, verbose=True)

    print(">>> Solving using Newton <<<")
    # ...
    PDE = PDE_newton
    if PDE.Dirichlet:
        U = PDE.unknown_dirichlet
    else:
        U = PDE.unknown
    # ...
    PDE_newton.solve(F, dF, u0=None, maxiter=100, rtol=1.e-6, verbose=True)

    print("norm using Picard  ", PDE_picard.norm(exact=u_exact))
    print("norm using Newton  ", PDE_newton.norm(exact=u_exact))
    # ...

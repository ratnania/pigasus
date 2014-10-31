# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from scipy.sparse.linalg import spsolve
from .basicPDE import *

__all__ = ['basicPDE_picard']

class basicPDE_picard(basicPDE):
    """
    A multidimentional nonlinear Poisson class solver using Picard algorithm.
        >>> import caid.cad_geometry  as cg
        >>> from caid.cad_geometry import line
        >>> import pylab                as pl


    """

    #: Doc comment for class attribute gallery.basicPDE.
    #: It can have multiple lines.
    def __init__(self, *args, **kwargs):
        """Creates a nonlinear basicPDE PDE solver based on Picard algorithm.

        geometry:
            The geometry must be an object cad_geometry.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.gallery.basicPDE`.

        """

        # ...

        basicPDE.__init__(self, *args, **kwargs)
        self.forceAssembly  = False
        self.Assembled      = False
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
        """
        Assemblies the linear system to be solved.
        """
        if f is not None:
            self.F_V.set_func(f)

        self.forceAssembly = update

        if self.forceAssembly or not self.Assembled:
            basicPDE.assembly(self)
            self.Assembled = True
    #-----------------------------------

    #-----------------------------------
    def solve(self, F, u0=None, maxiter=100, rtol=1.e-6 \
              , verbose=False, update=False):
        """
        solves the nonlinear basicPDE equation using PIcard algorithm
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
        self.assembly(update=update)

        # project u0 onto the discrete vectorial space
        self.initialize(u0=u0)

        # ...
        PDE  = self
        un   = PDE.unknown
        rhs  = self.rhs
        if self.Dirichlet:
            U = self.unknown_dirichlet
        else:
            U = self.unknown
        # ...

        # ...
        def rhs_func(x,y):
            return F(U,x,y)
        rhs.set_func(rhs_func)
        # ...

        # ...
        list_Err = [1.e6]
        un_values = un.get()
        i = 0
        while (list_Err[-1] > rtol) and (i < maxiter):
            U_old_values = U.get()

            # assembly the right hand side
            rhs.reset()
            self.update()
            # solve and update unew
            basicPDE.solve(self, rhs)

            U_values = U.get()
            err = np.linalg.norm(U_values-U_old_values)
            list_Err.append(err)
            if verbose:
                print(i, ": ","   |F(x)| = ", list_Err[-1])
            i += 1
        list_Err = np.asarray(list_Err[1:])
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
        return [4. * exp (_U)]

    def dF (U,x, y):
        U = U.evaluate()
        return[-4 * exp (U)]
    # ...

    # ...
    tc = {}

    tc['A'] = lambda x,y : [1., 0., 0., 1.]
    tc['u'] = lambda x,y : [0.]
    tc['f'] = lambda x,y : [0.]
    tc['AllDirichlet'] = True
    # ...

    PDE_picard = basicPDE_picard(  geometry=geo \
                                 , testcase=tc )

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

    print("norm using Picard  ", PDE_picard.norm(exact=u_exact))
    # ...

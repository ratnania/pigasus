# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from scipy.sparse.linalg import spsolve
from pigasus.fem.basicPDE import *

__all__ = ['bilaplacian']

class bilaplacian(basicPDE):
    """
    """

    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, geometry):
        """Creates an poisson PDE solver. arguments are the same as pigasus.__init__

        geometry:
            The geometry must be an object cad_geometry.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.basicPDE`.

        """

        # ...
        dim = geometry.dim
        if dim == 1:
            func_one  = lambda x : [ 1. ]
            func_zero = lambda x : [ 0. ]
            func_bip  = lambda x : [ 1. ]
        if dim == 2:
            func_one  = lambda x,y : [ 1. ]
            func_zero = lambda x,y : [ 0. ]
            func_bip  = lambda x,y : [ 1., 0., 0. \
                                     , 0., 0., 0. \
                                     , 0., 0., 1. ]
        if dim == 3:
            func_one  = lambda x,y,z : [ 1. ]
            func_zero = lambda x,y,z : [ 0. ]
            func_bip  = lambda x,y,z : [ 1., 0., 0. \
                                       , 0., 1., 0. \
                                       , 0., 0., 1. ]
        # ...

        testcase = {}

        testcase['D2'] = func_bip
        testcase['u']  = func_zero
        testcase['f']  = func_one

        testcase['AllDirichlet']  = True

        basicPDE.__init__(self, geometry=geometry, testcase=testcase)

        # ...
    #-----------------------------------

    #-----------------------------------
    def assembly(self, f=None, u=None):
        """
        Assemblies the linear system to be solved.
        """
        if f is not None:
            self.F_V.set_func(f)

        if u is not None:
            self.N_U.set_func(u)

        basicPDE.assembly(self)
    #-----------------------------------


if __name__ == '__main__':
    import caid.cad_geometry  as cg
    from caid.cad_geometry import line, bilinear
    import matplotlib.pyplot    as plt
    import numpy                as np
    from __main__ import __file__ as filename

    # ...
    sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt ; pi = np.pi
    # ...

    def testcase_square_Dirichlet():
        kx = 2. * pi ; ky = 2. * pi

        # exact solution
        u = lambda x,y : [sin ( kx * x ) * sin ( ky * y )]

        # rhs
        f = lambda x,y : [( kx**4 + ky**4 ) * sin ( kx * x ) * sin ( ky * y )]

        return f, u
    #-----------------------------------

    #-----------------------------------
    nx      = 31
    ny      = 31
    px      = 3
    py      = 3
    #-----------------------------------

    #-----------------------------------
    # ...
    from caid.cad_geometry import square as domain
    geo = domain(n=[nx,ny], p=[px,py])

    f, u = testcase_square_Dirichlet()
    PDE = bilaplacian(geometry=geo)
    # ...

    # ...
    PDE.assembly(f=f, u=u)
    PDE.solve()
    PDE.plot()  ; plt.colorbar()
    plt.savefig(filename.split('.py')[0]+'.png', format='png')

    normU = PDE.norm()
    print("norm U   = ", normU)
    # ...

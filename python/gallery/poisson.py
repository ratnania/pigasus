# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from pigasus.fem.basicPDE import *

__all__ = ['poisson']

class poisson(basicPDE):
    """
    A multidimentional EllipticPDE class solver.
        >>> import caid.cad_geometry  as cg
        >>> from caid.cad_geometry import line
        >>> import pylab                as pl
        >>> import numpy                as np

    creation of the geometry
        >>> nrb = line(p0=(0,0), p1=(1,0))
        >>> geo = cg.cad_geometry(geo=nrb)
        >>> geo.refine(id=0,list_p=[px-1])
        >>> tx = np.linspace(0.,1.,nx+2)[1:-1]
        >>> geo.refine(id=0, list_t=[tx])

    creation of the testcase: it must be a dictionary
        >>> testcase = {}

    we need to specify the list of faces for which we apply Homogeneous Dirichlet boundary condition
        >>> bc['list_DirFaces'] = [[1,2]]

    definition of the exact solution (if known) and the right hand side  (source/load term).
        >>> kx = 2. * pi
        >>> # exact solution
        >>> u = lambda x : sin ( kx * x )
        >>> # rhs
        >>> f = lambda x : ( kx**2) * sin ( kx * x )

    specifying the values of the solution on the boundary.
        >>> # values of u at the boundary
        >>> g1 = lambda x: u(x)
        >>> g2 = lambda x: u(x)
        >>> bc['list_faces'] = [[1,2]]
        >>> bc['list_g'] = [[g1, g2]]

    let's call pigasus now!
        >>> PDE = poisson(geometry=geo, bc=bc)
        >>> PDE.assembly(f=f, u=u)
        >>> PDE.solve()
        >>> PDE.plot(); pl.show()
        >>> PDE.norm()
        >>> PDE.get()

    """

    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, geometry, bc_dirichlet=None, bc_neumann=None,
                 AllDirichlet=None, Dirichlet=None, metric=None, solverInfo=None):
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

        testcase = {}

        testcase['A'] = func_stiff
#        testcase['b'] = func_one
        testcase['u'] = func_zero
        testcase['f'] = func_one
#        testcase['b'] = lambda x,y : [1.e-4 * 1.]

        if bc_dirichlet is not None:
            testcase['bc_dirichlet'] = bc_dirichlet

        if bc_neumann is not None:
            testcase['bc_neumann'] = bc_neumann

        if AllDirichlet is not None:
            testcase['AllDirichlet']  = AllDirichlet

        if Dirichlet is not None:
            testcase['Dirichlet']  = Dirichlet

        if metric is not None:
            testcase['metric']  = metric

        if solverInfo is not None:
            testcase['solverInfo']  = solverInfo

        basicPDE.__init__(self, geometry=geometry, testcase=testcase)

        self.forceAssembly  = False
        self.Assembled      = False
        # ...
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
    #-----------------------------------


if __name__ == '__main__':
    import caid.cad_geometry  as cg
    from caid.cad_geometry import line, bilinear
    import pylab                as pl
    import numpy                as np

    # ...
    sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt ; pi = np.pi
    # ...

    def testcase_line():
        kx = 2. * pi
        # exact solution
        u = lambda x : [sin ( kx * x )]
        # ...

        # rhs
        f = lambda x : [( kx**2) * sin ( kx * x )]
        # ...
        return f, u

    def testcase_square_Dirichlet():
        kx = 2. * pi ; ky = 2. * pi

        # exact solution
        u = lambda x,y : [sin ( kx * x ) * sin ( ky * y )]

        # rhs
        f = lambda x,y : [( kx**2 + ky**2 ) * sin ( kx * x ) * sin ( ky * y )]

        return f, u
    #-----------------------------------

    #-----------------------------------
    niter 	= 5000
    nx      = 15
    ny      = 15
    px      = 2
    py      = 2
    #-----------------------------------

    #-----------------------------------
    # ...
    from caid.cad_geometry import line
    geo1 = line(n=[nx], p=[px])

    f1, u1 = testcase_line()
    PDE1 = poisson(geometry=geo1)
    # ...

    # ...
    from caid.cad_geometry import square as domain
    geo2 = domain(n=[nx,ny], p=[px,py])

    f2, u2 = testcase_square_Dirichlet()
    PDE2 = poisson(geometry=geo2)
    # ...

    # ...
    PDE1.assembly(f=f1, u=u1)
    PDE2.assembly(f=f2, u=u2)
    # ...

    # ...
    PDE1.solve()
    PDE2.solve()
    # ...

    # ...
    PDE1.plot()  ; pl.show()
    PDE2.plot()  ; pl.show()
    # ...

    # ...
    normU1 = PDE1.norm()
    print("norm U-1D   = ", normU1)
    # ...

    # ...
    normU2 = PDE2.norm()
    print("norm U-2D   = ", normU2)
    # ...

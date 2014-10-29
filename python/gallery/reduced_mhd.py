# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from pigasus.fem.basicPDE import multi_basicPDE

__all__ = ['reduced_mhd']

class reduced_mhd(multi_basicPDE):
    """
    A multidimentional nonlinear Poisson class solver using Picard algorithm.
        >>> import caid.cad_geometry  as cg
        >>> from caid.cad_geometry import line
        >>> import pylab                as pl


    """

    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, *args, **kwargs):
        """Creates a one-step parabolic PDE solver.

        geometry:
            The geometry must be an object cad_geometry.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.gallery.poisson`.

        """

        # ...
        parabolic.__init__(self, *args, **kwargs)
        # ...

    def solve(self, niter):
        """
        updates directly the field of the implicit Operator
        """
        E = self.list_PDE[0]
        I = self.list_PDE[1]

        # ...
        un   = E.unknown
        unew = I.unknown

#        un.set(E.rhs)
        # ...

        # ...
        for i in range(0,niter):

            rhs = E.dot(un)
            I.solve(rhs)

            un.set(unew)
        # ...
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
#        # ...
#        bc = {}
#
#        # ...
#        # values of u at the boundary
#        # ...
#        g1 = lambda x: 0.
#        g2 = lambda x: 0.
#
#        bc['list_faces'] = [[1,2]]
#        bc['list_g'] = [[g1, g2]]
#        # ...

        kx = 2. * pi
        u = lambda x : [0.] # will be set for the implicit part
        # i -> implicit part, e -> explicit part
        # constant term
        fe = lambda x : [( kx**2) * sin ( kx * x )]
        fi = lambda x : [0.]
        # mass
        be = lambda x : [1.]
        bi = lambda x : [1.]
        # advection
        #ve = lambda x : [0., 0.]
        #vi = lambda x : [0., 0.]
        # transpose advection
        #tve = lambda x : [0., 0.]
        #tvi = lambda x : [0., 0.]
        # stiffness
        alpha = 1. ; dt = 1.e-3
        Ae = lambda x : [-dt*(1.-alpha)*1.]
        Ai = lambda x : [dt*alpha*1.]
        # dictionaries for implicit and explicit parts
        tci = {}       ; tce = {}
        tci['f'] = fi  ; tce['f'] = fe
        tci['b'] = bi  ; tce['b'] = be
        #tci['v'] = vi  ; tce['v'] = ve
        #tci['tv'] = tvi ; tce['tv'] = tve
        tci['A'] = Ai  ; tce['A'] = Ae
        tci['u'] = u
        tci['AllDirichlet'] = True ; tce['AllDirichlet'] = True

        return tce, tci

    #-----------------------------------
    niter 	= 500
    nx      = 31
    ny      = 31
    px      = 2
    py      = 2
    #-----------------------------------

    #-----------------------------------
    # ...
    from caid.cad_geometry import line
    geo1 = line(n=[nx], p=[px])

    tce1, tci1 = testcase_line()
    PDE1 = multi_basicPDE(geo1, tce1, tci1)
    # ...

    # ...
    PDE1.assembly()
    # ...

    # ...
    PDE1.solve(niter)
    # ...

    # ...
    PDE1.plot()  ; pl.show()
    # ...
#
#    # ...
#    normU1 = PDE1.norm()
#    print "norm U-1D   = ", normU1
#    # ...
#

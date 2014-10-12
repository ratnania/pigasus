# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from pigasus.fem.basicPDE import *

__all__ = ['parabolic', 'onestep']

class parabolic(object):
    """
    A multidimentional Parabolic equations class solver.
        >>> import igakit.cad_geometry  as cg

    """

    #: Doc comment for class attribute gallery.poisson.
    #: It can have multiple lines.
    def __init__(self, geometry, list_tc, bc_dirichlet=None, bc_neumann=None,
                 AllDirichlet=None, Dirichlet=None, metric=None, solverInfo=None):
        """Creates an poisson PDE solver. arguments are the same as pigasus.__init__

        geometry:
            The geometry must be an object cad_geometry.

        list_tc:
            contains the testcases list, needed for each implicit/explicit operator. The user must first specify the explicit terms, then the implicit ones.

        Returns:
           A PDE object.

        .. note::

            See also: :ref:`fem.basicPDE`.

        """

        # ...
        self.dim = geometry.dim
        self.nPDEs = len(list_tc)
        # ...

        # ...
        self.list_PDE = []
        for tc in list_tc:
            testcase = {}

            try:
                testcase['A'] = tc['A']
            except:
                pass

            try:
                testcase['v'] = tc['v']
            except:
                pass

            try:
                testcase['w'] = tc['w']
            except:
                pass

            try:
                testcase['b'] = tc['b']
            except:
                pass

            try:
                testcase['D2'] = tc['D2']
            except:
                pass

            try:
                testcase['u'] = tc['u']
            except:
                pass

            try:
                testcase['f'] = tc['f']
            except:
                pass

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

            # ...
            PDE = basicPDE(geometry=geometry, testcase=testcase)
            self.list_PDE.append(PDE)
            # ...
        # ...
    #-----------------------------------

    #-----------------------------------
    def __del__(self):
        for PDE in self.list_PDE:
            PDE.__del__()
    #-----------------------------------

    #-----------------------------------
    def free(self):
        for PDE in self.list_PDE:
            PDE.free()
    #-----------------------------------

    #-----------------------------------
    def initialize(self, list_u0=None):
        _list_u0        = [None]*self.nPDEs
        if list_u0 is not None:
            _list_u0    = list_u0

        for u0,PDE in zip(_list_u0,self.list_PDE):
            U = PDE.unknown
            if u0 is None:
                U.set(np.zeros(U.size))
            else:
                PDE.interpolate(u0, field=U)
    #-----------------------------------

    #-----------------------------------
    def assembly(self, list_f=None, list_update=None):
        _list_f         = [None]*self.nPDEs
        if list_f is not None:
            _list_f     = list_f

        _list_update    = [False]*self.nPDEs
        if list_update is not None:
            _list_update = list_update

        for f,update,PDE in zip(_list_f, _list_update, self.list_PDE):
            if f is not None:
                self.F_V.set_func(f)

            PDE.forceAssembly = update

            if PDE.forceAssembly or not PDE.Assembled:
                PDE.assembly()
    #-----------------------------------

    #-----------------------------------
    def plot(self):
        PDE = self.list_PDE[-1]
        PDE.plot()
    #-----------------------------------

class onestep(parabolic):
    """
    A multidimentional nonlinear Poisson class solver using Picard algorithm.
        >>> import igakit.cad_geometry  as cg
        >>> from igakit.cad import line
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
    import igakit.cad_geometry  as cg
    from igakit.cad import line, bilinear
    import pylab                as pl
    import numpy                as np

    # ...
    sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt ; pi = np.pi
    # ...

    def testcase_line():
        # ...
        bc = {}

        # ...
        # values of u at the boundary
        # ...
        g1 = lambda x: 0.
        g2 = lambda x: 0.

        bc['list_faces'] = [[1,2]]
        bc['list_g'] = [[g1, g2]]
        # ...

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
        dicti = {}       ; dicte = {}
        dicti['f'] = fi  ; dicte['f'] = fe
        dicti['b'] = bi  ; dicte['b'] = be
        #dicti['v'] = vi  ; dicte['v'] = ve
        #dicti['tv'] = tvi ; dicte['tv'] = tve
        dicti['A'] = Ai  ; dicte['A'] = Ae
        dicti['u'] = u

        return bc, dicte, dicti

    #-----------------------------------
    niter 	= 500
    nx      = 31
    ny      = 31
    px      = 2
    py      = 2
    #-----------------------------------

    #-----------------------------------
    # ...
    from igakit.cad_geometry import line
    geo1 = line(n=[nx], p=[px])

    bc1, dicte1, dicti1 = testcase_line()
    PDE1 = parabolic(geo1, dicte1, dicti1, bc=bc1)
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

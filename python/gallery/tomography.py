# -*- coding: UTF-8 -*-
#! /usr/bin/python

import matplotlib.pyplot    as plt
import numpy                as np
from time import time
from pigasus.fem.basicPDE import basicPDE
from caid.cad_geometry import cad_geometry, cad_nurbs
from caid.cad_geometry import square as patch
from caid.core.bspline import bsp
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
from pigasus.fit.utils import *

#-----------------------------------
class component(object):

    def __init__(self):
        # ... global indices of basis related to this component
        self._basis_indices = []

    @property
    def basis_indices(self):
        return self._basis_indices
#-----------------------------------

#-----------------------------------
class component_geometry_limiter(component):

    def __init__(self, geometry):
        self._geometry = geometry
        component.__init__(self)

    @property
    def basis_indices(self):
        return self._basis_indices

    @property
    def geometry(self):
        return self._geometry
#-----------------------------------

#-----------------------------------
class tomography(object):
    def __init__(self, geometry, mu=None, alpha=1., rational=0):
        """
        """
        self._geometry              = geometry
        self._postAssembly          = False
        self._alpha                 = alpha
        self._mu                    = mu
        self._rational              = rational
        self._list_PDE              = []
        self._entropy               = None
        self._total_power           = None
        self._fisher_information    = None

        # ...
        dim = geometry.dim
        if dim == 1:
            func_mass   = lambda x : [ 1. ]
            func_stiff  = lambda x : [ 1. ]
            func_bip    = lambda x : [ 1. ]
        if dim == 2:
            func_mass   = lambda x,y : [  1. ]
            func_stiff  = lambda x,y : [  1., 0. \
                                        , 0., 1. ]
            func_bip    = lambda x,y : [ 1., 0., 0. \
                                     , 0., 2., 0. \
                                     , 0., 0., 1. ]
        if dim == 3:
            # TODO
            func_mass   = lambda x,y,z : [  1. ]
            func_stiff  = lambda x,y,z : [ 1., 0., 0. \
                                         , 0., 1., 0. \
                                         , 0., 0., 1. ]
            func_bip    = lambda x,y,z : [ None ]
        # ...

        # ... reference PDE
        testcase = {}
        PDE = basicPDE(geometry=geometry, testcase=testcase)
        self._list_PDE.append(PDE)
        PDE_ref = PDE
        self._PDE_ref = PDE
        # ...

        # ... PDE for the mass operator
        testcase = {}
        testcase['b'] = func_mass
        PDE = basicPDE(geometry=geometry, testcase=testcase, V=PDE_ref.V)
        self._list_PDE.append(PDE)
        # ...

        # ... PDE for the stiffness operator
        testcase = {}
        testcase['A'] = func_stiff
        PDE = basicPDE(geometry=geometry, testcase=testcase, V=PDE_ref.V)
        self._list_PDE.append(PDE)
        # ...

        # ... PDE for the bilaplacian operator
        testcase = {}
        testcase['D2'] = func_bip
        PDE = basicPDE(geometry=geometry, testcase=testcase, V=PDE_ref.V)
        self._list_PDE.append(PDE)
        # ...

        # ... assembly all PDE
        for PDE in self._list_PDE:
            PDE.assembly()
        # ...

    @property
    def space(self):
        return self._PDE_ref.space

    @property
    def ID_loc(self):
        return self._PDE_ref.space.connectivity.ID_loc

    @property
    def entropy(self):
        return self._entropy

    @property
    def total_power(self):
        return self._total_power

    @property
    def fisher_information(self):
        return self._fisher_information
#-----------------------------------


if __name__ == "__main__":
    from caid.cad_geometry import square as domain
    nx = 15 ; ny = 15
    px = 2  ; py = 2
    geo = domain(n=[nx,ny], p=[px,py])

    ToMo = tomography(geo)

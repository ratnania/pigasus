# -*- coding: UTF-8 -*-
#! /usr/bin/python

import sys
import numpy as np
from scipy.sparse.linalg import spsolve

from pigasus.fem.constants import *
from pigasus.fem.field import *
from pigasus.fem.norm import *
from pigasus.fem.grids import *
from pigasus.fem.matrix import *
from pigasus.fem.space import *
from pigasus.fem.blockdata import BlockMatrix, BlockVector
from pigasus.fem.pigasus import *

__all__ = ['EllipticPDEMG']

class EllipticPDEMG(pigasus):
    def __init__(self, *args, **kwargs):
        """
        constructs the object EllipticPDE with
        fem                     : the fem main object for pigasus
        geo                     : the cad_geometry object
        testcase                : the current test case, from Monge-Ampere solutions
        """
        pigasus.__init__(self, *args, **kwargs)

        # ...
        try:
            self.PDE_H = self.testcase['PDE_H']
        except:
            print("PDE_H must be specified")
            sys.exit(0)
        # ...

        # ...
        try:
            self.PDE_h = self.testcase['PDE_h']
        except:
            print("PDE_h must be specified")
            sys.exit(0)
        # ...

        #-----------------------------------
        self.geometry_H = self.PDE_H.geometry
        self.geometry_h = self.PDE_h.geometry

        nrb = self.geometry_h[0]
        dim = self.geometry_h.dim
        list_n = nrb.shape
        list_p = nrb.degree

        lpi_ordregl = list_p
        #-----------------------------------

        #-----------------------------------
        # space declaration :coarse
        #-----------------------------------
        V_H = space(geometry=self.geometry_H)
        V_H.dirichlet(faces=self.testcase['list_DirFaces'])
        V_H.set_boundary_conditions()
        V_H.create_grids(type="legendre", k=lpi_ordregl)
        #-----------------------------------

        #-----------------------------------
        # space declaration : fine
        #-----------------------------------
        V_h = space(geometry=self.geometry_h)
        V_h.dirichlet(faces=self.testcase['list_DirFaces'])
        V_h.set_boundary_conditions()
        V_h.create_grids(type="legendre", k=lpi_ordregl)
        #-----------------------------------

        #-----------------------------------
        # fields
        #-----------------------------------
        if dim == 1:
            func_one   = lambda x : [ 1. ]
            func_zero  = lambda x : [ 0. ]
        if dim == 2:
            func_one   = lambda x,y : [ 1. ]
            func_zero  = lambda x,y : [ 0. ]
        if dim == 3:
            func_one   = lambda x,y,z : [ 1. ]
            func_zero  = lambda x,y,z : [ 0. ]
        # ...

        F_H = field(space=V_H, func = func_one)
        F_h = field(space=V_h, func = func_one)
        #-----------------------------------

        #-----------------------------------
        # Differential Operators
        #-----------------------------------
        #-----------------------------------

        #-----------------------------------
        # Save access for data
        #-----------------------------------
        self.V_H    = V_H
        self.V_h    = V_h
        self.F_H    = F_H
        self.F_h    = F_h
        #-----------------------------------

        self.initInterpolation()
        self.initRestriction()
    #-----------------------------------

    #-----------------------------------
    def assembly(self):
        pass
    #-----------------------------------

    #-----------------------------------
    def initInterpolation(self):
        from scipy.sparse import csr_matrix
        dim = self.geometry_h.dim
        if dim ==1:
            from .splineRefMat import constructCurveMatrix as constructMatrix
        if dim ==2:
            from .splineRefMat import constructSurfaceMatrix as constructMatrix
        if dim ==3:
            print("initInterpolation: Not yet implemented for 3D")

        nrb_H = self.geometry_H[0]
        nrb_h = self.geometry_h[0]

        if dim ==1:
            knots_H = nrb_H.knots[0]
            knots_h = nrb_h.knots[0]

            n = nrb_H.shape[0]
            p = nrb_H.degree[0]

            list_r = [r for r in knots_h if r not in knots_H]

            M = constructMatrix(list_r, p, n, knots_H, bc=0)
            self.PH_h = M

#                self.PH_h = csr_matrix(M.todense()[1:-1,1:-1])

        if dim ==2:
            u_H1,u_H2   = nrb_H.knots
            n_H1,n_H2   = nrb_H.shape
            p_H1,p_H2   = nrb_H.degree

            u_h1,u_h2   = nrb_h.knots

            list_r1 = [r for r in u_h1 if r not in u_H1]
            list_r2 = [r for r in u_h2 if r not in u_H2]

            M, [n,m] = constructMatrix(  list_r1, list_r2 \
                                , p_H1, p_H2 \
                                , n_H1, n_H2 \
                                , u_H1, u_H2 \
                                , bc=[0,0])
            self.PH_h = M

#        print 'M.shape : ', M.shape
#        print 'self.PH_h.shape : ', self.PH_h.shape
    #-----------------------------------

    #-----------------------------------
    def initRestriction(self):
        list_H, list_h = self.compute_H()
        r = list_h[0]/list_H[0]
        self.Ph_H = self.PH_h.transpose().tocsr()
        self.Ph_H *= r
#        print 'self.Ph_H.shape : ', self.Ph_H.shape
    #-----------------------------------

    #-----------------------------------
    def interpolation(self, vH):
        vh = self.PH_h.dot(vH)
        return vh
    #-----------------------------------

    #-----------------------------------
    def restriction(self, vh):
        vH = self.Ph_H.dot(vh)
        return vH
    #-----------------------------------

    #-----------------------------------
    def compute_H(self):
        geo_H = self.geometry_H
        geo_h = self.geometry_h

        dim = geo_H.dim

        list_H = []
        list_h = []

        for (nrb_H, nrb_h) in zip(geo_H, geo_h):
            H = 1.
            h = 1.
            for d in range(0,dim):
                p_H = nrb_H.degree[d]
                p_h = nrb_H.degree[d]
                H *= nrb_H.knots[d][p_H+1]-nrb_H.knots[d][0]
                h *= nrb_h.knots[d][p_h+1]-nrb_h.knots[d][0]

            list_H.append(H)
            list_h.append(h)

        return list_H, list_h
    #-----------------------------------

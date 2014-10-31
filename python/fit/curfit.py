# -*- coding: UTF-8 -*-
#! /usr/bin/python

import matplotlib.pyplot    as plt
import numpy                as np
from time import time
from pigasus.fem.basicPDE import basicPDE
from caid.cad_geometry import cad_geometry, cad_nurbs
from caid.cad_geometry import square as patch
from caid.core.bspline import bsp
#from igakit.nurbs import NURBS
from scipy.io import mmwrite
from scipy.sparse import coo_matrix
from pigasus.fit.utils import *

#-----------------------------------
class curfit(object):
    def __init__(self, geometry, uk=None, PDE=None, constraints=[],
                 mu=None, alpha=1., rational=0):
        """
        initialize the curfit object
        PDE is the Differential operator to use for smoothing (usually a 2nd
        order)
        constraints is a list of dictionaries that must be of the following form
        constraints[i] is  {'patch_id_m', 'face_m', 'patch_id_s', 'face_s',
        'type'}
        patch_id_m is the master patch id
        face_m     is the face id in the master patch
        patch_id_s is the slave  patch id
        face_s     is the face id in the slave  patch
        type       is the constraint's type: C1, C2, ... (default: C1)
        ib         is the starting index in the face element (default:0 )
        ie         is the ending   index in the face element (default:-1 )
        """
        self.geometry       = geometry
        self.postAssembly   = False
        self.nConstraints   = 0
        self.ConstIndices   = []
        self.ConstVals      = []
        self.ConstRHSs      = []
        self.constraints    = constraints
        self.alpha          = alpha
        self.mu             = mu
        self.rational       = rational
        self.Rd             = 2 # TODO make it automatic
        self.RHS            = None
        self.Mat            = None

        if PDE is not None:
            self.PDE = PDE
        else:
            from pigasus.fem.basicPDE import basicPDE

            func_zero = lambda x : [ 0.]
            func_bip  = lambda x : [ 1.]
            testcase = {}

            testcase['D2'] = func_bip
            testcase['u']  = func_zero
            testcase['f']  = func_zero

            self.PDE = basicPDE(geometry=geometry, testcase=testcase)

        # assembly the PDE
        self.PDE.assembly()

        self.ID_loc  = self.PDE.space.connectivity.ID_loc
        if uk is not None:
            self.Dt, self.Mat = self.updateMatrix(uk)
            self.postAssembly = True

    @property
    def system(self):
        return self.PDE.system.get()

    @property
    def space(self):
        return self.PDE.space

    # ...
    def updateMatrix(self, lists_uk):
        geo = self.geometry
        ID_loc = self.ID_loc
        Mat = self.system

        Mat_shape = Mat.shape
    #   print "shape ", Mat_shape
        for patch_id in range(0, geo.npatchs):
            nrb = geo[patch_id]
            list_uk = lists_uk[patch_id]

            ID = ID_loc[patch_id]
            # ...

#            tb = time()
            list_basis = evalBasis1D(nrb, list_uk, rational=self.rational)
#            print "basis    ",list_basis
#            te = time()
#            print "elapsed time basis evaluation ", te-tb

#            tb = time()
            Dt, Mk = addContributions1D(nrb, list_uk, list_basis, ID, Mat_shape)
#            te = time()
#            print "elapsed time discrete Mass ", te-tb

            if self.mu is None:
                Mat_norm = np.linalg.norm(Mat.todense())
                Mk_norm  = np.linalg.norm(Mk.todense())
                mu = Mk_norm / Mat_norm
            else:
                mu = self.mu

            alpha = self.alpha

            Mat = Mk + alpha * mu * Mat

        return Dt, Mat
    # ...

    def fit(self, xyzk, uk=None):
        if (not self.postAssembly) and (uk is None):
            print("You must run updateMatrix before!")
            return

        if (not self.postAssembly) and (uk is not None):
            self.Dt, self.Mat = self.updateMatrix(uk)
            self.updateGlobalSystem()
            self.postAssembly = True

        from scipy.sparse.linalg import cg as solve
        N_ini   = self.system.shape[0]
        N_final = N_ini+self.nConstraints
        cxyz = []
        for i,xk in enumerate(xyzk):
            x  = np.zeros(N_final)
            x[:N_ini] = self.Dt * xk
            # ... update with rhs for some Constraints
            x += np.asarray(self.RHS[i])
            # ...
#            np.savetxt("x"+str(id(xk))+".txt", x)
            Cx = solve(self.Mat, x)[0]
            cx = Cx[:N_ini]
            cxyz.append(cx)
#        from scipy.linalg import det
#        print "det = ", det(self.Mat.todense())
#        print "---"
#        print ">> xk "
#        print xk
#        print "---"
#        print ">> self.Dt "
#        print self.Dt
#        print "---"
#        print ">> x "
#        print x
#        print "---"
#        print self.Mat
#        print "---"
#        print Cx
#        print "---"
        return cxyz

    def construct(self, xyzk, uk=None):
        cxyz = self.fit(xyzk, uk=uk)

        U = self.PDE.unknown

        geo_ini = self.geometry

        geo_f = cad_geometry()
        for patch_id in range(0, self.geometry.npatchs):
            nrb = self.geometry[patch_id]

            C = np.zeros(list(nrb.shape)+[3])
            for i,cx in enumerate(cxyz):
                U.set(cx)
                C[...,i] =  U.tomatrix(patch_id).reshape(nrb.shape)

            srf = cad_nurbs(nrb.knots, C, weights=nrb.weights)
            srf.orientation = nrb.orientation
            srf.rational = nrb.rational
            geo_f.append(srf)

        geo_f._internal_faces = geo_ini._internal_faces
        geo_f._external_faces = geo_ini._external_faces
        geo_f._connectivity   = geo_ini._connectivity

        return geo_f

    # ...
    def genLineIndices(self, patch_id, face, shift=0):
        if face == 0:
            list_A = np.asarray([self.ID_loc[patch_id][0+shift]])
        if face == 1:
            list_A = np.asarray([self.ID_loc[patch_id][-1-shift]])

        list_A = list_A - 1 # 0 based index
        return list_A
    # ...

    # ...
    def addConstraint(self, patch_id_m, fm, typeC \
                      , patch_id_s=None, fs=None \
                      , list_vals=None):

        # ... C0 Condition
        if typeC == "C0":
            nrb_m = self.geometry[patch_id_m] # master
            list_Am = self.genLineIndices(patch_id_m, fm, shift=0)

            list_As = None
            if (patch_id_s is not None) and (fs is not None):
                nrb_s = self.geometry[patch_id_s] # slave
                list_As = self.genLineIndices(patch_id_s, fs, shift=0)

            rows, values, rhs = genLineC0Constraint( self.Rd \
                                                   , list_Am \
                                                   , list_As=list_As \
                                                   , list_vals=list_vals \
                                                   )
        # ...

        # ... C1 Condition
        if typeC == "C1":
            nrb_m = self.geometry[patch_id_m] # master
            list_Am = self.genLineIndices(patch_id_m, fm, shift=0)
            list_Bm = self.genLineIndices(patch_id_m, fm, shift=1)
            cm      = computeC1Coef1D(nrb_m, fm)

            list_As = None
            list_Bs = None
            cs      = None
            if (patch_id_s is not None) and (fs is not None):
                nrb_s = self.geometry[patch_id_s] # slave

                list_As = self.genLineIndices(patch_id_s, fs, shift=0)
                list_Bs = self.genLineIndices(patch_id_s, fs, shift=1)
                cs      = computeC1Coef1D(nrb_s, fs)

            rows, values, rhs = genLineC1Constraint( self.Rd \
                                                   , list_Am \
                                                   , list_Bm \
                                                   , cm \
                                                   , list_As=list_As \
                                                   , list_Bs=list_Bs \
                                                   , cs=cs \
                                                   , list_vals=list_vals \
                                                   )
        # ...

        self.nConstraints += len(rows)
        return rows, values, rhs
    # ...

    # ...
    def updateGlobalSystem(self):
        # ...
        for constraint in self.constraints:
            patch_id_m  = None
            fm          = None
            typeC       = "C1"
            patch_id_s  = None
            fs          = None
            list_vals   = []
#            print "Constraint ", constraint
            for key, value in constraint.items():
#                print key
#                print value

                if key == "patch_id_m":
                    patch_id_m  = int(value)

                if key == "face_m":
                    fm  = int(value)

                if key == "type":
                    typeC  = str(value)

                if key == "patch_id_s":
                    patch_id_s  = int(value)

                if key == "face_s":
                    fs  = int(value)

                if key == "values":
                    list_vals = value

#            print "==============="
#            print "patch_id_s ", patch_id_s
#            print "fs ", fs
#            print "list_vals ", list_vals
#            print "==============="

            if len(list_vals) == 0:
                list_vals = None

            rows, values, rhs = self.addConstraint( patch_id_m, fm \
                                                  , typeC \
                                                  , patch_id_s=patch_id_s, fs=fs \
                                                  , list_vals=list_vals \
                                                 )

            self.ConstIndices   += rows
            self.ConstVals      += values
#            print "r  ", rhs
            self.ConstRHSs.append(rhs)
        # ...

        # ...
        Mat = self.Mat.tocoo()
        n,m = Mat.shape
#        print "Mat.shape ", self.Mat.shape

        allData = list(Mat.data)
        allRows = list(Mat.row)
        allCols = list(Mat.col)

        allRHS = []
        for i in range(0,self.Rd):
            allRHS.append(list(np.zeros(m)))
#        print "allRHS ",  allRHS
#        print "self.ConstRHSs  ",  self.ConstRHSs
#        print "len ", len(allRHS[0])

        nCurrent = n
        for indices, data, rhs in zip(  self.ConstIndices \
                                      , self.ConstVals \
                                      , self.ConstRHSs):
            allData += data
            allRows += indices
            allCols += [nCurrent]*len(indices)
            # treatment of the transposed part
            allData += data
            allCols += indices
            allRows += [nCurrent]*len(indices)
            # treatment of the RHS
#            print "rhs ", rhs
            _rhs = [list(L) for L in zip(*rhs)]
#            print ">>> _rhs : ",_rhs
#            print "len ", len(allRHS[0])
            for i in range(0, self.Rd):
#                print ">> i ", i
#                print "_rhs[i] ", _rhs[i]
                allRHS[i] += _rhs[i]
#            print "len ", len(allRHS[0])
            nCurrent += 1

        n += self.nConstraints
        m += self.nConstraints
        shp = [n,m]

        self.Mat = coo_matrix((allData, (allRows, allCols)), shape=shp)
        self.Mat.tocsr()
#        print "allRHS ",  allRHS
#        print "allRHS-x ",  allRHS[0]
#        print "allRHS-y ",  allRHS[1]
#        print len(allRHS[0])
#        print len(allRHS[1])
        self.RHS = [np.asarray(rhs) for rhs in allRHS]
#        print "*********"
#        print "Mat.shape ", self.Mat.shape
#        print "RHS.shape ", [R.shape for R in self.RHS]
#        print "*********"
#        from scipy.io import mmwrite
#        mmwrite("Mat.mtx", self.Mat)
#        np.savetxt("RHS.txt", self.RHS)
        # ...

    # ...
#-----------------------------------

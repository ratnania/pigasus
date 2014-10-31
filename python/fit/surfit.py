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
class surfit(object):
    def __init__(self, geometry, uvk=None, PDE=None, constraints=[],
                 mu=None, alpha=1., rational=0):
        """
        initialize the surfit object
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

            func_zero = lambda x,y : [ 0. ]
            func_bip  = lambda x,y : [ 1., 0., 0. \
                                     , 0., 2., 0. \
                                     , 0., 0., 1. ]
            testcase = {}

            testcase['D2'] = func_bip
            testcase['u']  = func_zero
            testcase['f']  = func_zero

            self.PDE = basicPDE(geometry=geometry, testcase=testcase)

        # assembly the PDE
        self.PDE.assembly()

        self.ID_loc  = self.PDE.space.connectivity.ID_loc
        if uvk is not None:
            self.Dt, self.Mat = self.updateMatrix(uvk)
            self.postAssembly = True

    @property
    def system(self):
        return self.PDE.system.get()

    @property
    def space(self):
        return self.PDE.space

    # ...
    def updateMatrix(self, lists_uvk, verbose=False):
        lists_uk = lists_uvk[0]
        lists_vk = lists_uvk[1]
        geo = self.geometry
        ID_loc = self.ID_loc
        Mat = self.system

        Mat_shape = Mat.shape
#        print "shape ", Mat_shape

        # ... construct the D matrix for the first patch
        patch_id = 0
        nrb = geo[patch_id]
        list_uk = lists_uk[patch_id]
        list_vk = lists_vk[patch_id]
        ID = ID_loc[patch_id]
#        np.savetxt("ID"+str(patch_id)+".txt",np.asarray(ID), fmt='%d')

        list_basis = evalBasis2D(nrb, list_uk, list_vk \
                                 , rational=self.rational \
                                 , verbose=verbose)
        rows, cols, data = addContributions2D(nrb, list_uk, list_vk \
                                    , list_basis, ID \
                                    , verbose=verbose)
#        np.savetxt("rows"+str(patch_id)+".txt", rows, fmt='%d')
#        np.savetxt("cols"+str(patch_id)+".txt", cols, fmt='%d')
        # ...

        # ... coompute and update the D matrix for the other patchs
        for patch_id in range(1, geo.npatchs):
            nrb = geo[patch_id]
            list_uk = lists_uk[patch_id]
            list_vk = lists_vk[patch_id]

            ID = ID_loc[patch_id]
#            np.savetxt("ID"+str(patch_id)+".txt",np.asarray(ID), fmt='%d')

            # ...

            list_basis = evalBasis2D(nrb, list_uk, list_vk \
                                     , rational=self.rational \
                                     , verbose=verbose)
            rowsk, colsk, datak = addContributions2D(nrb, list_uk, list_vk \
                                        , list_basis, ID \
                                        , verbose=verbose)

            colsk = np.asarray(colsk)
            colsk += np.max(cols)+1
#            np.savetxt("rows"+str(patch_id)+".txt", rowsk, fmt='%d')
#            np.savetxt("cols"+str(patch_id)+".txt", colsk, fmt='%d')
            rows = np.concatenate([rows, rowsk])
            cols = np.concatenate([cols, colsk])
            data = np.concatenate([data, datak])

        ntotal = len(np.concatenate(lists_uk))
        shp = [Mat_shape[0],ntotal]

        D   = coo_matrix((data, (rows, cols)), shape=shp)
        D   = D.tocsr()
        Dt  = D.transpose().tocsr()
        DtD = D * Dt
        Mk = DtD.transpose().tocsr()
#        print "Mat_shape    ", Mat_shape
#        print "system shape ", self.system.shape
#        print "Mk shape     ", Mk.shape

        if self.mu is None:
            Mat_norm = np.linalg.norm(Mat.todense())
            Mk_norm  = np.linalg.norm(Mk.todense())
            mu = Mk_norm / Mat_norm
        else:
            mu = self.mu

        alpha = self.alpha

        Mat = Mk + alpha * mu * Mat
#        from scipy.io import mmwrite
#        mmwrite("Mat.mtx", Mat)
#        mmwrite("D.mtx", D)

        return D, Mat
    # ...

    def fit(self, xyzk, uvk=None):
        if (not self.postAssembly) and (uvk is None):
            print("You must run updateMatrix before!")
            return

        if (not self.postAssembly) and (uvk is not None):
            self.Dt, self.Mat = self.updateMatrix(uvk)
            self.postAssembly = True

        self.updateGlobalSystem()

        from scipy.sparse.linalg import cg as solve
        N_ini   = self.system.shape[0]
        N_final = N_ini+self.nConstraints
#        print "N_ini ", N_ini
#        print "N_final ", N_final
#        print "nconstraints ", self.nConstraints
        cxyz = []
        # xk is a list of arrays. one array per patch
        for i,xk in enumerate(xyzk):
#            print "==== i "+str(i)+" ==="
            x  = np.zeros(N_final)
            # we must transform xk into one single array
            _xk = np.concatenate(xk)
#            print "xk ", len(_xk)
#            print "Dt ", self.Dt.shape
            x[:N_ini] = self.Dt * _xk
            # ... update with rhs for some Constraints
            x += np.asarray(self.RHS[i])
            # ...
#            np.savetxt("x"+str(id(xk))+".txt", x)
            Cx = solve(self.Mat, x)[0]
            cx = Cx[:N_ini]
#            np.savetxt("xk"+str(i)+".txt", _xk)
#            np.savetxt("x"+str(i)+".txt", x)
#            np.savetxt("xa"+str(i)+".txt", x[:81])
#            np.savetxt("xb"+str(i)+".txt", x[81:])
#            np.savetxt("cx"+str(i)+".txt", cx)
#            print ">> x "
#            print x
#            print "---"
#            print "---"
#            print cx
#            print "---"
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
#        print self.Mat
        return cxyz


    def construct(self, xyzk, uvk=None, exportGeometry=True):
        cxyz = self.fit(xyzk, uvk=uvk)

        U = self.PDE.unknown

        geo_ini = self.geometry
        if exportGeometry:
            geo_f = cad_geometry()
        else:
            list_nrb = []

        for patch_id in range(0, self.geometry.npatchs):
            nrb = self.geometry[patch_id]

            C = np.zeros(list(nrb.shape)+[3])
            for i,cx in enumerate(cxyz):
                U.set(cx)
                C[...,i] =  U.tomatrix(patch_id).reshape(nrb.shape)

            if exportGeometry:
                srf = cad_nurbs(nrb.knots, C, weights=nrb.weights)
                srf.orientation = nrb.orientation
                srf.rational = nrb.rational
                geo_f.append(srf)
            else:
                srf = cad_nurbs(nrb.knots, C, weights=nrb.weights)
                list_nrb.append(srf)

        if exportGeometry:
            geo_f.set_internal_faces(geo_ini.internal_faces)
            geo_f.set_external_faces(geo_ini.external_faces)
            geo_f.set_connectivity(geo_ini.connectivity)
            return geo_f
        else:
            return list_nrb

    # ...
    def genLineIndices(self, patch_id, face, shift=0, ib=0, ie=None):
        if face == 0:
            if ie is None:
                list_A = self.ID_loc[patch_id][ib:,0+shift]
            else:
                list_A = self.ID_loc[patch_id][ib:ie,0+shift]
        if face == 1:
            if ie is None:
                list_A = self.ID_loc[patch_id][0+shift,ib:]
            else:
                list_A = self.ID_loc[patch_id][0+shift,ib:ie]
        if face == 2:
            if ie is None:
                list_A = self.ID_loc[patch_id][ib:,-1-shift]
            else:
                list_A = self.ID_loc[patch_id][ib:ie,-1-shift]
        if face == 3:
            if ie is None:
                list_A = self.ID_loc[patch_id][-1-shift,ib:]
            else:
                list_A = self.ID_loc[patch_id][-1-shift,ib:ie]

        list_A = list_A - 1 # 0 based index
        return list_A
    # ...

    # ...
    def addConstraint(self, patch_id_m, fm, typeC \
                      , patch_id_s=None, fs=None \
                      , list_vals=None \
                      , ib=0, ie=None):

        # ... C0 Condition
        if typeC == "C0":
            nrb_m = self.geometry[patch_id_m] # master
            list_Am = self.genLineIndices(patch_id_m, fm, shift=0)
#            print "list_Am ", list_Am

            list_As = None
            if (patch_id_s is not None) and (fs is not None):
                nrb_s = self.geometry[patch_id_s] # slave
                list_As = self.genLineIndices(patch_id_s, fs, shift=0)
#                print "list_As ", list_As

            rows, values, rhs = genLineC0Constraint( self.Rd \
                                                   , list_Am \
                                                   , list_As=list_As \
                                                   , list_vals=list_vals \
                                                   )
#            print "rows   ", rows
#            print "values ", values
#            print "rhs    ", rhs
        # ...

        # ... C1 Condition
        if typeC == "C1":
            nrb_m = self.geometry[patch_id_m] # master
            list_Am = self.genLineIndices(patch_id_m, fm, shift=0, ib=ib, ie=ie)
            list_Bm = self.genLineIndices(patch_id_m, fm, shift=1, ib=ib, ie=ie)
            cm      = computeC1Coef2D(nrb_m, fm)

            list_As = None
            list_Bs = None
            cs      = None
            if (patch_id_s is not None) and (fs is not None):
                nrb_s = self.geometry[patch_id_s] # slave

                list_As = self.genLineIndices(patch_id_s, fs, shift=0, ib=ib, ie=ie)
                list_Bs = self.genLineIndices(patch_id_s, fs, shift=1, ib=ib, ie=ie)
                cs = computeC1Coef2D(nrb_s, fs)

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
            ib          = 0
            ie          = -1
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

                if key == "ib":
                    ib  = int(value)

                if key == "ie":
                    ie  = int(value)

#            print "==============="
#            print "patch_id_m ", patch_id_m
#            print "patch_id_s ", patch_id_s
#            print "fm ", fm
#            print "fs ", fs
#            print "list_vals ", list_vals
#            print "==============="

            if len(list_vals) == 0:
                list_vals = None

            rows, values, rhs = self.addConstraint(  patch_id_m, fm \
                                              , typeC \
                                              , patch_id_s=patch_id_s, fs=fs \
                                              , list_vals=list_vals \
                                              , ib=ib \
                                              , ie=ie)

            self.ConstIndices   += rows
            self.ConstVals      += values
            self.ConstRHSs      += rhs
#            print "r  ", rhs
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

#        print "self.ConstIndices  ",  self.ConstIndices
#        print "self.ConstVals  ",  self.ConstVals
#        print "self.ConstRHSs  ",  self.ConstRHSs

        nCurrent = n
        for indices, data, rhs in zip(  self.ConstIndices \
                                      , self.ConstVals \
                                      , self.ConstRHSs):
#            print "indices ", indices
#            print "data    ", data
#            print "rhs     ", rhs
            allData += data
            allRows += indices
#            allCols += nCurrent*np.ones(len(indices), dtype=np.int) #[nCurrent]*len(indices)
            allCols += [nCurrent]*len(indices)
            # treatment of the transposed part
            allData += data
            allCols += indices
#            allRows += nCurrent*np.ones(len(indices), dtype=np.int) #[nCurrent]*len(indices)
            allRows += [nCurrent]*len(indices)
            # treatment of the RHS
#            print "rhs ", rhs
            _rhs = rhs #[list(L) for L in zip(*rhs)]
#            print ">>> _rhs : ",_rhs
#            print "len ", len(allRHS[0])
            for i in range(0, self.Rd):
#                print ">> i ", i
#                print "_rhs[i] ", _rhs[i]
                allRHS[i].append(_rhs[i])
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
#        print "Final Mat.shape ", self.Mat.shape
#        print "Final RHS.shape ", [R.shape for R in self.RHS]
#        print "*********"
#        from scipy.io import mmwrite
#        mmwrite("Mat.mtx", self.Mat)
#        np.savetxt("RHS.txt", self.RHS)
        # ...

    # ...
#-----------------------------------

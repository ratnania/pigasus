# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['fast_iga']
__date__ ="$Feb 11, 2012 11:51:45 PM$"

from scipy.sparse.linalg import splu
import numpy as np
from .pigasusObject import *

class fast_iga(pigasusObject):
    def __init__(self, matrices=[[]], useSPLU=True):
        pigasusObject.__init__(self)

        self.matrices = matrices
        self.nmatrices = len(matrices)
        self.com.usefiga = True

        # setting the number of matrices
        self.com.pyfem.set_nmatrices_figa(self.nmatrices)
        self.com.pyfem.create_figa_partone()

        # setting left and right matrices
        LEFT  = 1
        RIGHT = 0
        li_i = 0
        for [M,K] in matrices:
            li_i += 1
            self.com.pyfem.set_matrix_figa(M.id, li_i, LEFT)
            self.com.pyfem.set_matrix_figa(K.id, li_i, RIGHT)

        self.com.pyfem.initialize_figa()

        lpi_n = self.com.pyfem.get_nmnk_figa()
#        print "lpi_n=", lpi_n
        self.nM = lpi_n[0]
        self.nK = lpi_n[1]

        # if we use super-LU we need to import all matrices, and factorize them
        # first we import all matrices in self.list_Matrix
        self.list_Matrix = []
        for li_i in range(0, self.nK):
            lo_M = self.export_matrix(li_i+1)
            self.list_Matrix.append(lo_M)
        # factorization
        self.list_opM = []
        for M in self.list_Matrix:
            self.list_opM.append(splu(M.tocsc()))

    def to_common_basis(self):
        self.com.pyfem.fwdfft_figa()

    def from_common_basis(self):
        self.com.pyfem.bwdfft_figa()

    def export_matrix (self, ai_id):
        """
        this routine exports temp matrices that are used in the Fast-IGA, in CSR format
        """
        import numpy as np
        import scipy.sparse as ss

        li_nR = 0
        li_nC = 0
        li_nel = 0
        li_nR, li_nC, li_nel = self.com.pyfem.getparamcsrmatrix_figa ( ai_id )
        lpr_a  = np.zeros(li_nel)
        lpi_ia = np.zeros(li_nR+1)
        lpi_ja = np.zeros(li_nel)
        lpr_a, lpi_ia, lpi_ja = self.com.pyfem.getarrayparamcsrmatrix_figa ( ai_id, li_nR, li_nC, li_nel )

        #indices must begin from 0, and not 1 like in fortran
        lpi_ja = lpi_ja - 1
        lpi_ia = lpi_ia - 1
        return ss.csr_matrix((lpr_a, lpi_ja, lpi_ia) )

    def set_rhs(self, apr_F):
        self.com.pyfem.set_rhs_figa(apr_F)

    def get_rhs(self):
        return self.com.pyfem.get_rhs_figa(self.nM, self.nK)

    def set_u(self, apr_U):
        self.com.pyfem.set_sol_figa(apr_U)

    def get_u(self):
        return self.com.pyfem.get_sol_figa(self.nM, self.nK)

    def solve(self, apr_F):
        """
        *******************************************
                SOLVING THE LINEAR SYSTEM
                we solve each column-system
        *******************************************
        TODO attention a la parite
        ... destinguish between the cases n even or odd
        """
        li_n = self.nM
        lpr_F1D = np.zeros(li_n, dtype=np.double)
        lpr_U1D = np.zeros(li_n, dtype=np.double)
        apr_U   = np.zeros((self.nK,li_n), dtype=np.double)

        li_mmax = (self.nK / 2) - 1
        #... m is even
        if self.nK % 2 :
            #... treatment of the 0-mode
            lpr_F1D = apr_F[0, :]
            lpr_U1D = self.list_opM[0].solve ( lpr_F1D )
            apr_U[0, :] = lpr_U1D

            for li_j in range(1,li_mmax+1):
                #... treatment of the mode 2j-1
                lpr_F1D = apr_F[2 * li_j - 1, :]
                lpr_U1D = self.list_opM[li_j].solve ( lpr_F1D )
                apr_U[2 * li_j - 1, :] = lpr_U1D
                #... treatment of the mode 2j
                lpr_F1D = apr_F[2 * li_j, :]
                lpr_U1D = self.list_opM[li_j].solve ( lpr_F1D )
                apr_U[2 * li_j, :] = lpr_U1D

            #... treatment of the last mode
            lpr_F1D = apr_F[self.nK-1, :]
            lpr_U1D = self.list_opM[li_mmax+1].solve ( lpr_F1D )
            apr_U[self.nK-1, :] = lpr_U1D
        #... m is odd
        else :
             #... treatment of the 0-mode
            lpr_F1D = apr_F[0, :]
            lpr_U1D = self.list_opM[0].solve ( lpr_F1D )
            apr_U[0, :] = lpr_U1D

            for li_j in range(1,( ( self.nK - 1 ) / 2 ) + 1):
                #... treatment of the mode 2j-1
                lpr_F1D = apr_F[2 * li_j - 1, :]
                lpr_U1D = self.list_opM[li_j].solve ( lpr_F1D )
                apr_U[2 * li_j - 1, :] = lpr_U1D
                #... treatment of the mode 2j
                lpr_F1D = apr_F[2 * li_j, :]
                lpr_U1D = self.list_opM[li_j].solve ( lpr_F1D )
                apr_U[2 * li_j, :] = lpr_U1D

        return apr_U

    def decomptensor(self, V, W, X, ai_patch, apr_values):
        li_dimV = V.geometry[ai_patch].dim
        li_dimW = W.geometry[ai_patch].dim
        li_dimX = X.geometry[ai_patch].dim
        if (li_dimV==3) and (li_dimW==1) and (li_dimX==2) :
            return self.decomptensor_312(V, W, X, ai_patch, apr_values)
        if (li_dimV==2) and (li_dimW==1) and (li_dimX==1) :
            return self.decomptensor_211(V, W, X, ai_patch, apr_values)

    def comptensor(self, V, W, X, ai_patch, apr_values):
        li_dimV = V.geometry[ai_patch].dim
        li_dimW = W.geometry[ai_patch].dim
        li_dimX = X.geometry[ai_patch].dim
        if (li_dimV==3) and (li_dimW==1) and (li_dimX==2) :
            return self.comptensor_312(V, W, X, ai_patch, apr_values)
        if (li_dimV==2) and (li_dimW==1) and (li_dimX==1) :
            return self.comptensor_211(V, W, X, ai_patch, apr_values)

    def decomptensor_312(self, V, W, X, ai_patch, apr_values):
        """
        this routine returns the input vector apr_values into a matrix
        apr_values is a vector on the space V
        the output is a matrix on (W,X)
        dim V = 3
        dim W = 1
        dim X = 2
        """
        lpi_nV = V.geometry[ai_patch].shape
        lpi_nW = W.geometry[ai_patch].shape
        lpi_nX = X.geometry[ai_patch].shape

        import numpy as np
        lpr_out = np.zeros((W.connectivity.size,X.connectivity.size),dtype=np.double)

        if (lpi_nV[0]==lpi_nW[0]) and (lpi_nV[1]==lpi_nX[0]) and (lpi_nV[2]==lpi_nX[1]) :
            li_ni = lpi_nV[0]
            li_nj = lpi_nV[1]
            li_nk = lpi_nV[2]
            for li_k in range(0, li_nk):
                for li_j in range(0, li_nj):
                    for li_i in range(0, li_ni):
                        li_A_V = V.geometry[ai_patch].get_A_ind([li_i,li_j,li_k])
                        li_A_W = W.geometry[ai_patch].get_A_ind([li_i])
                        li_A_X = X.geometry[ai_patch].get_A_ind([li_j,li_k])
                        li_P_V = V.connectivity.ID[li_A_V]
                        li_P_W = W.connectivity.ID[li_A_W]
                        li_P_X = X.connectivity.ID[li_A_X]
                        if (li_P_V!=0) and (li_P_W!=0) and (li_P_X!=0):
                            lpr_out[li_P_W-1,li_P_X-1] = apr_values[li_P_V-1]
        else:
            print("Error : you must verify spaces dimensions")
            print("lpi_nV =", lpi_nV)
            print("lpi_nW =", lpi_nW)
            print("lpi_nX =", lpi_nX)

        return lpr_out

    def comptensor_312(self, V, W, X, ai_patch, apr_values):
        """
        this routine returns the input vector apr_values into a matrix
        apr_values is a matrix on (W,X)
        the output is a vector on the space V
        dim V = 3
        dim W = 1
        dim X = 2
        """
        lpi_nV = V.geometry[ai_patch].shape
        lpi_nW = W.geometry[ai_patch].shape
        lpi_nX = X.geometry[ai_patch].shape

        import numpy as np
        lpr_out = np.zeros(V.connectivity.size,dtype=np.double)

        if (lpi_nV[0]==lpi_nW[0]) and (lpi_nV[1]==lpi_nX[0]) and (lpi_nV[2]==lpi_nX[1]) :
            li_ni = lpi_nV[0]
            li_nj = lpi_nV[1]
            li_nk = lpi_nV[2]
            for li_k in range(0, li_nk):
                for li_j in range(0, li_nj):
                    for li_i in range(0, li_ni):
                        li_A_V = V.geometry[ai_patch].get_A_ind([li_i,li_j,li_k])
                        li_A_W = W.geometry[ai_patch].get_A_ind([li_i])
                        li_A_X = X.geometry[ai_patch].get_A_ind([li_j,li_k])
                        li_P_V = V.connectivity.ID[li_A_V]
                        li_P_W = W.connectivity.ID[li_A_W]
                        li_P_X = X.connectivity.ID[li_A_X]
                        if (li_P_V!=0) and (li_P_W!=0) and (li_P_X!=0):
                            lpr_out[li_P_V-1] = apr_values[li_P_W-1,li_P_X-1]
        else:
            print("Error : you must verify spaces dimensions")
            print("lpi_nV =", lpi_nV)
            print("lpi_nW =", lpi_nW)
            print("lpi_nX =", lpi_nX)

        return lpr_out

    def decomptensor_211(self, V, W, X, ai_patch, apr_values):
        """
        this routine returns the input vector apr_values into a matrix
        apr_values is a vector on the space V
        the output is a matrix on (W,X)
        dim V = 2
        dim W = 1
        dim X = 1
        """
        lpi_nV = V.geometry[ai_patch].shape
        lpi_nW = W.geometry[ai_patch].shape
        lpi_nX = X.geometry[ai_patch].shape

        import numpy as np
        lpr_out = np.zeros((W.connectivity.size,X.connectivity.size),dtype=np.double)

        if (lpi_nV[0]==lpi_nW[0]) and (lpi_nV[1]==lpi_nX[0]) :
            li_ni = lpi_nV[0]
            li_nj = lpi_nV[1]
            for li_j in range(0, li_nj):
                for li_i in range(0, li_ni):
                    li_A_V = V.geometry[ai_patch].get_A_ind([li_i,li_j])
                    li_A_W = W.geometry[ai_patch].get_A_ind([li_i])
                    li_A_X = X.geometry[ai_patch].get_A_ind([li_j])
#                    print "************"
#                    print "li_i,li_j=", li_i,li_j
#                    print "li_A_V=", li_A_V
#                    print "li_A_W,li_A_X=", li_A_W,li_A_X
                    li_P_V = V.connectivity.ID[li_A_V]
                    li_P_W = W.connectivity.ID[li_A_W]
                    li_P_X = X.connectivity.ID[li_A_X]
#                    print "li_P_V=", li_P_V
#                    print "li_P_W,li_P_X=", li_P_W,li_P_X
                    if (li_P_V!=0) and (li_P_W!=0) and (li_P_X!=0):
                        lpr_out[li_P_W-1,li_P_X-1] = apr_values[li_P_V-1]
        else:
            print("Error : you must verify spaces dimensions")
            print("lpi_nV =", lpi_nV)
            print("lpi_nW =", lpi_nW)
            print("lpi_nX =", lpi_nX)
        return lpr_out

    def comptensor_211(self, V, W, X, ai_patch, apr_values):
        """
        this routine returns the input vector apr_values into a matrix
        apr_values is a matrix on (W,X)
        the output is a vector on the space V
        dim V = 2
        dim W = 1
        dim X = 1
        """
        lpi_nV = V.geometry[ai_patch].shape
        lpi_nW = W.geometry[ai_patch].shape
        lpi_nX = X.geometry[ai_patch].shape

        import numpy as np
        lpr_out = np.zeros(V.connectivity.size,dtype=np.double)

        if (lpi_nV[0]==lpi_nW[0]) and (lpi_nV[1]==lpi_nX[0]) :
            li_ni = lpi_nV[0]
            li_nj = lpi_nV[1]
            for li_j in range(0, li_nj):
                for li_i in range(0, li_ni):
                    li_A_V = V.geometry[ai_patch].get_A_ind([li_i,li_j])
                    li_A_W = W.geometry[ai_patch].get_A_ind([li_i])
                    li_A_X = X.geometry[ai_patch].get_A_ind([li_j])
#                    print "************"
#                    print "li_i,li_j=", li_i,li_j
#                    print "li_A_V=", li_A_V
#                    print "li_A_W,li_A_X=", li_A_W,li_A_X
                    li_P_V = V.connectivity.ID[li_A_V]
                    li_P_W = W.connectivity.ID[li_A_W]
                    li_P_X = X.connectivity.ID[li_A_X]
#                    print "li_P_V=", li_P_V
#                    print "li_P_W,li_P_X=", li_P_W,li_P_X
                    if (li_P_V!=0) and (li_P_W!=0) and (li_P_X!=0):
                        lpr_out[li_P_V-1] = apr_values[li_P_W-1,li_P_X-1]
        else:
            print("Error : you must verify spaces dimensions")
            print("lpi_nV =", lpi_nV)
            print("lpi_nW =", lpi_nW)
            print("lpi_nX =", lpi_nX)
        return lpr_out






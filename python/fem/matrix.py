# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['matrix']
__date__ ="$Jan 11, 2012 3:32:21 PM$"

from . pigasusObject import *
#from graph import *
from . import common_obj as _com
from . import constants as _cst
from numpy import zeros
##############################################################################
#
#       matrix Class
#
##############################################################################
class matrix(pigasusObject):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
#        obj = object.__new__(typ, *args, **kwargs)
        obj.id = None
        return obj

    def __init__ ( self, graph=None ):
        pigasusObject.__init__(self)

        self.operators  = []
        self.IJVAssembly = False

        self.graph = graph
        self.type = _cst.GENERIC_MATRIX
        if graph is None:
            self.type = _cst.BLOCK_MATRIX

        self.id = self.com.nmatrices
        self.com.nmatrices += 1
        self.com.matrices.append(self)

        self.setInfoData()

    def setInfoData(self):
        """
        prints informations about the current matrix
        """
        try:
            self.infoData['id'] = str(self.id)
        except:
            self.infoData['id'] = None

        try:
            self.infoData['operators'] = str([O.id for O in self.operators])
        except:
            self.infoData['operators'] = None

#        try:
#            self.infoData['shape'] = str(self.shape)
#        except:
#            self.infoData['shape'] = None

    def append_operator(self, O):
        self.operators.append(O)

    def get(self):
        """
        returns a scipy.sparse.csr_matrix
        """
        import numpy as np
        if self.IJVAssembly:
            li_size = 0
            li_size = self.com.pyfem.pyfem_getsizeijv ( self.id )
            lpi_i = np.zeros(li_size, dtype=np.int)
            lpi_j = np.zeros(li_size, dtype=np.int)
            lpr_v = np.zeros(li_size)
            lpi_i, lpi_j, lpr_v = self.com.pyfem.pyfem_getdataijv ( self.id, li_size )
            return lpi_i, lpi_j, lpr_v

        else:
            import scipy.sparse as ss

            li_nR = 0
            li_nC = 0
            li_nel = 0
            li_nR, li_nC, li_nel = self.com.pyfem.getparamcsrmatrix ( self.id )
            lpr_a  = np.zeros(li_nel)
            lpi_ia = np.zeros(li_nR+1)
            lpi_ja = np.zeros(li_nel)
            lpr_a, lpi_ia, lpi_ja = self.com.pyfem.getarrayparamcsrmatrix ( self.id, li_nR, li_nel )

            return ss.csr_matrix((lpr_a, lpi_ja, lpi_ia) )

    def set(self, other):
        """
        sets the current matrix-coeff to a given value or array
        """

        # ...
        def _createfromscipy(A):
            """
            creates a matrix from a csr description
            A is a scipy.sparse matrix
            """
            [li_nR, li_nC] = A.get_shape()
            li_nnz = A.getnnz()
            self.com.pyfem.createcsrmatrix( self.id, li_nC  \
            , A.indptr, A.indices, li_nR, li_nnz)
        # ...

        # ...
        def _setfromscipy(A):
            """
            creates a matrix from a csr description
            A is a scipy.sparse matrix
            """
            li_nnz = A.getnnz()
            self.com.pyfem.setcsrmatrix( self.id, A.data, li_nnz)
        # ...

        # ...
        if self.id is None:
            self._values = other
        else:
#            if isFloat(other):
#                self._set_val(other)
#            if isNumpyArray(other):
#                self._set_array(other)
            if _com.isMatrix(other):
                self._set_matrix(other)
            if _com.isScipyMatrix(other):
                other = other.tocsr()
                _createfromscipy(other)
                _setfromscipy(other)


    def reset(self):
        self.com.pyfem.pyfem_reset_matrix ( self.id )

    def assembly(self):
        """
        assembling the current matrix will run assemble all its operators
        """
        # TODO
        pass

    def solve(self, rhs):
        """
        solves the linear system Ax=y
        A:
            the current matrix
        x,y: may be fields or numpy arrays
        """
        # TODO
        pass

    @property
    def nnz(self):
        li_nR = 0
        li_nC = 0
        li_nel = 0
        li_nR, li_nC, li_nel = self.com.pyfem.getparamcsrmatrix ( self.id )
        return li_nel

    @property
    def shape(self):
        li_nR = 0
        li_nC = 0
        li_nel = 0
        li_nR, li_nC, li_nel = self.com.pyfem.getparamcsrmatrix ( self.id )
        return [li_nR, li_nC]

#    def clone(self):
#        M = matrix.__new__(matrix)
##        M.id = None
#        M.set(self.get ())
#        return M
#
#    def zeros_like(self):
#        M = matrix.__new__(matrix)
##        M.id = None
#        M.set(0.0 * self.get ())
#        return M

    def save(self, filename):
        from scipy.io import mmwrite
        mmwrite(filename, self.get())

    def load(self, filename):
        from scipy.io import mmread
        self.set(mmread(filename).tocsr())

    def dot(self, other):
        """
        scalar dot between a matrix and a field
        """
        return self.__mul__(other)

#    def __radd__(self, other):
#        return self.__add__(other)
#
    def __rmul__(self, other):
        return self.__mul__(other)

    def __pos__(self):
        self.__mul__(1.0)

    def __neg__(self):
        self.__mul__(-1.0)

#    def __sub__(self, other):
#        return self.__add__(-other)
#
#    def __isub__(self, other):
#        return self.__iadd__(-other)
#
#    def __iadd__(self, other):
#        if _com.isFloat(other):
#            self.com.pyfem.matrix_add_scal (self.id, other, self.id)
#        if _com.isMatrix(other):
#            self.com.pyfem.matrix_add_matrix (self.id,other.id,self.id)
#        if _com.isScipyMatrix(other):
#            self.set(other + self.get ())
#        return self
#
    def __imul__(self, other):
        if _com.isFloat(other):
            self.com.pyfem.matrix_mult_scal (self.id, other, self.id)
        if _com.isMatrix(other):
            self.com.pyfem.matrix_mult_matrix (self.id,other.id,self.id)
        if _com.isScipyMatrix(other):
            self.set(other * self.get ())
        return self

#    def __add__(self, other):
#        if _com.isFloat(other):
#            M = matrix.__new__(matrix)
#            M.set(other + self.get ())
#        elif _com.isMatrix(other):
#            M = matrix.__new__(matrix)
#            M.set(other.get () + self.get ())
#        elif _com.isScipyMatrix(other):
#            M = matrix.__new__(matrix)
#            M.set(other + self.get ())
#        return M
#
#        print "Not yet implemented in __add__"

    def __mul__(self, other):
        if _com.isNumpyArray(other):
#            print "mult matrix with numpy array"
            n,m = self.shape
            x = zeros(n)
            x = self.com.pyfem.matrix_mult_array2(self.id, other, n )
            return x
        if _com.isField(other):
            from . import field as fi
            F = fi.field.__new__(fi.field)
            n,m = self.shape
            F.set(self.com.pyfem.matrix_mult_array(self.id, other.id, m ))
            return F
        if _com.isFloat(other):
            M = matrix.__new__(matrix)
            M.set(other * self.get ())
            return M
        if _com.isMatrix(other):
            M = matrix.__new__(matrix)
            M.set(other.get () * self.get ())
            return M
        if _com.isScipyMatrix(other):
            M = matrix.__new__(matrix)
            M.set(other * self.get ())
            return M

        print("Not yet implemented in __mul__")

##############################################################################

if __name__ == '__main__':
#    import fem      as fem
#    fe = fem.fem()

    print("test Matrix class")
    M = matrix()

#    fe.initialize()
#    print M
#    print "---"
#
#    from scipy import sparse
#    from numpy import array
#    I = array([0,3,1,2,3,3,2])
#    J = array([0,3,1,2,2,1,1])
#    V = array([4,5,7,9,-1,1,-2])
##    A = sparse.coo_matrix((V,(I,J)),shape=(4,4)).tocsr()
##
##    M.set(A)
#    print "done"

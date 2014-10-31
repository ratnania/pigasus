# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['oper']
__date__ ="$Jan 11, 2012 3:32:21 PM$"

from .pigasusObject import *
from . import common_obj as _com
from . import constants as _cst
import numpy as _np
##############################################################################
#
#       oper Class
#
##############################################################################
class oper(pigasusObject):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
#        obj = object.__new__(typ, *args, **kwargs)
        obj.id = None
        return obj

    def __init__ ( self, spaces, type, func = None \
    , parameval = False, paramevalfunc = False \
    , transpose = False, func_arguments = [], addto = [] ):
        """
        parameval :
            True if we want to have the eval-basis on the parametric domain : grad,...
        transpose :
            if we want to transpose the matrix during the assembling process
        """

        pigasusObject.__init__(self)

        self.func_arguments = func_arguments
        self.paramevalfunc  = paramevalfunc
        self.type           = type

        self.nparam         = 0
        self.parameval      = parameval
        self.transpose      = transpose
        self.spaces         = spaces
        self.list_addto     = []

        # master sp is the id of space which will have more informations
        # among all spaces related to the matrix
        # default = 0, unless one of the spaces is composed
        li_master_sp        = 0
        self.space          = spaces[li_master_sp]
        self.loc_id         = spaces[li_master_sp].grids.add_operator_id(self)
        self.dim            = spaces[li_master_sp].dim

        self._set_nparam()
        self._set_nderiv()

        if func is None:
            func = self._set_default_func()

        # ...
        self.set_func(func)
        # ...

        self.id = self.com.noperators
        self.com.noperators += 1
        self.com.operators.append(self)

        # TODO : mettre le setInfoDATA dans l'initialisation de fem
#        self.setInfoData()

    def set_func(self, func):
        """
        this sets the param-function of the current operator
        """
        from . import func_tools as ft
        if len(self.func_arguments) == 0 :
            from .utils import function
            self.func = function(func, space=self.space)
        else:
            raise("Not used anymore. Dead code")

    def setInfoData(self):
        """
        prints informations about the current matrix
        """
        try:
            self.infoData['id'] = str(self.id)
        except:
            self.infoData['id'] = None

        try:
            self.infoData['type'] = str(self.type)
        except:
            self.infoData['type'] = None

        try:
            self.infoData['nparam'] = str(self.nparam)
        except:
            self.infoData['nparam'] = None

        try:
            self.infoData['nderiv'] = str(self.nderiv)
        except:
            self.infoData['nderiv'] = None

        try:
            self.infoData['spaces'] = str([V.id for V in self.spaces])
        except:
            self.infoData['spaces'] = None

        try:
            txt = "\n"
            for data in self.list_addto:
                txt += "        matrix : " + str(data[0].id) + " scale : " + str(data[1]) + "\n"
            self.infoData['addto'] = txt
        except:
            self.infoData['addto'] = None

    def addto(self, M, scale=1.0):
        self.list_addto.append([M, scale])
        M.append_operator(self)

    def _set_default_func(self):
        """
        sets the number of parameters for the current matrix
        """
        dim = self.dim
        if ( self.type in [ _cst.MASS ] ):
            if dim == 1:
                func = lambda x : [ 1. ]
            if dim == 2:
                func = lambda x,y : [ 1. ]
            if dim == 3:
                func = lambda x,y,z : [ 1. ]
            return func
        if ( self.type in [ _cst.STIFFNESS ] ):
            if dim == 1:
                func = lambda x : [ 1. ]
            if dim == 2:
                func = lambda x,y : [  1., 0. \
                                     , 0., 1. ]
            if dim == 3:
                func = lambda x,y,z : [ 1., 0., 0. \
                                      , 0., 1., 0. \
                                      , 0., 0., 1. ]
            return func
        if ( self.type in [ _cst.ADVECTION ] ):
            if dim == 1:
                func = lambda x : [ 1. ]
            if dim == 2:
                func = lambda x,y : [ 1., 1. ]
            if dim == 3:
                func = lambda x,y,z : [ 1., 1., 1. ]
            return func
        if ( self.type in [ _cst.SECOND_DERIV ] ):
            if dim == 1:
                func = lambda x : [ 1. ]
            if dim == 2:
                func = lambda x,y,z : [ 1., 0., 0. \
                                      , 0., 0., 0. \
                                      , 0., 0., 1. ]
            if dim == 3:
                func = lambda x,y,z : [ 1., 1., 1., 0., 0., 0. ]
                raise("Not yet implemented")
            return func
        else :
            print("oper _set_default_func : type not implemented yet")
            import sys; sys.exit(1)

    def _set_nparam(self):
        """
        sets the number of parameters for the current matrix
        """

        if ( self.type in [ _cst.MASS ] ):
             self.nparam = 1
             return
        if ( self.type in [ _cst.STIFFNESS ] ):
            li_dim = self.space.dim
            self.nparam = li_dim**2
            return
        if ( self.type in [ _cst.ADVECTION ] ):
            li_dim = self.space.dim
            self.nparam = li_dim
            return
        if ( self.type in [ _cst.SECOND_DERIV ] ):
            li_dim = self.space.dim
            if li_dim == 1 :
                self.nparam = 1
            if li_dim == 2 :
                self.nparam = 9
            if li_dim not in [1,2] :
                print("Matrix : _set_nparam , attention popur le moment on gere que le 1D et 2D")
                raise("Not yet implemented")
                self.nparam = 3
            return
        if ( self.type in [ _cst.COMPOSITION ] ):
            self.nparam = 0
            return
        else :
            print("MATRIX-_set_nparam : type not implemented yet")
            import sys; sys.exit(1)

    def _set_nderiv(self):
        """
        sets the number of derivative involved in the current matrix
        """

        if ( self.type in [ _cst.MASS, _cst.ADVECTION, _cst.STIFFNESS ] ):
             self.nderiv = 1
             return

        if ( self.type in [ _cst.SECOND_DERIV ] ):
            self.nderiv = 2
            return
        if ( self.type in [ _cst.COMPOSITION ] ):
            self.nderiv = 0
            return

        self.nderiv = 0
        return

    def evalfunc(self, ai_patch, apr_points, elts=None):
        """
        Evaluation of the param-function over a given list of points
        """
        V = self.space

        lpr_pts = V.get_points()
        list_pts = []
        for i in range(0, V.dim):
            list_pts.append(lpr_pts[i,0,:])
        lpr_pts = list_pts

#        print "=============================="
#        print "%%% patch-id ", V.currentPatchID
#        print "%%% lpr_pts.shape ", len(lpr_pts)
#        print "%%% id-space ", id(V)
#        print "%%% id-field ", id(self)
        values  = self.func(lpr_pts)
#        print "%%% values.len ", len(values)
#        print "%%% values.shape ", values.shape
        return values
#        pvalues = self.pfunc(V.get_sites())
#        n,m = values.shape
#        v = _np.zeros((1,m))
#        for i in range(0,n):
#            v += values[i,:] * pvalues[i,:]
#
#        return v

#        # test if we are in the non linear case
#        if len(self.func_arguments) == 0:
##                            print "classic"
#            if not self.paramevalfunc :
#                lpr_val = self._evalfunc_std(ai_patch, apr_points, elts)
#            else:
#                lpr_parampts = self.space.get_parametricPoints(ai_patch_id=ai_patch)
#                lpr_val = self._evalfunc_std(ai_patch, lpr_parampts, elts)
##                            print "lpr_val.shape = ", lpr_val.shape
##                            print "lpr_val = ", lpr_val
#        else:
##            print "non classic"
#            lpr_val = self._evalfunc_generic(ai_patch,apr_points, elts)
##                        list_matrices[li_matrix].append(lpr_val)
##                            print "lpr_val = ", lpr_val
##                            print "lpr_val.shape = ", lpr_val.shape
##                        elapsed = (clock() - start)
##                        print ("CPU time for evaluating Fields and Matrices is :" + str (elapsed))
#        return lpr_val

    def _evalfunc_std(self, ai_patch, apr_points, elts=None):
        """
        sequential version of the evaluation
        """
        return self.func(apr_points)

    def varfunc(self, fields=[]):
        """
        this routine defines variables for the matrix function,
        variables must be of the form  (fields=[...])
        this is to implement functions of the form:
        F([u,v],[x,y]) where u and v are fields
        """
        self.func_arguments = fields

    def _evalfunc_generic(self, ai_patch, apr_points, elts=None):
        # using fields
        fields = self.func_arguments
        if len(fields) == 1 :
            F = fields[0]
#            print F.type, F.operator
            if ( (F.type == _cst.FIELD_OPERATOR) and (F.operator in [_cst.GRAD, _cst.CURL])):
#                print "related field : ", F.id
                lpr_val = F.eval(ai_patch)
            else :
#                print "Matrix : _evalfunc_generic for NON LINEAR evaluation"
                lpr_val = self._evalfunc_nonlin(ai_patch, apr_points, elts)
        else :
            lpr_val = self._evalfunc_nonlin(ai_patch, apr_points, elts)

        return lpr_val

    def _evalfunc_nonlin(self, ai_patch, apr_points, elts=None):
        """
        sequential version of the evaluation
        """
        # loop over fields involved in the function
        list_values_F = []
        for F in self.func_arguments:
            list_val = F.eval(ai_patch, elts)

            list_values_F.append(list_val)

        # TODO to change when passed to ndof > 1
        lpr_val = self.func(list_values_F, apr_points)
        return lpr_val

    def clone(self):
        # TODO
        pass

    def __radd__(self, other):
        print("TODO")
        pass
#        return self.__add__(other)

    def __rmul__(self, other):
        print("TODO")
        pass
#        return self.__mul__(other)

    def __pos__(self):
        print("TODO")
        pass
#        self.__mul__(1.0)

    def __neg__(self):
        print("TODO")
        pass
#        self.__mul__(-1.0)

    def __sub__(self, other):
        print("TODO")
        pass
#        return self.__add__(-other)

    def __isub__(self, other):
        print("TODO")
        pass
#        return self.__iadd__(-other)

    def __iadd__(self, other):
        print("TODO")
        pass
#        if _com.isFloat(other):
#            self.com.pyfem.matrix_add_scal (self.id, other, self.id)
#        if _com.isMatrix(other):
#            self.com.pyfem.matrix_add_matrix (self.id,other.id,self.id)
#        if _com.isScipyMatrix(other):
#            self.set(other + self.get ())
#        return self

    def __imul__(self, other):
        print("TODO")
        pass
#        if _com.isFloat(other):
#            self.com.pyfem.matrix_mult_scal (self.id, other, self.id)
#        if _com.isMatrix(other):
#            self.com.pyfem.matrix_mult_matrix (self.id,other.id,self.id)
#        if _com.isScipyMatrix(other):
#            self.set(other * self.get ())
#        return self

    def __add__(self, other):
        print("TODO")
        pass
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
        print("TODO")
        pass
#        if _com.isNumpyArray(other):
#            import field as fi
#            F = fi.field.__new__(fi.field)
#            F.set(self.com.pyfem.matrix_mult_array2(self.id, other, other.size ))
#            return F
#        if _com.isField(other):
#            import field as fi
#            F = fi.field.__new__(fi.field)
#            F.set(self.com.pyfem.matrix_mult_array(self.id, other.id, other.size ))
#            return F
#        if _com.isFloat(other):
#            M = matrix.__new__(matrix)
#            M.set(other * self.get ())
#            return M
#        if _com.isMatrix(other):
#            M = matrix.__new__(matrix)
#            M.set(other.get () + self.get ())
#            return M
#        if _com.isScipyMatrix(other):
#            M = matrix.__new__(matrix)
#            M.set(other * self.get ())
#            return M
#
#        print "Not yet implemented in __mul__"

##############################################################################

if __name__ == '__main__':
    from . import fem      as fem
    fe = fem.fem()

    print("done")

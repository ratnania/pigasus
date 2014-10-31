# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['norm']
__date__ ="$Feb 14, 2012 11:40:06 AM$"

from . import common_obj as _com
from . import constants as _cst
import numpy as _np
from .pigasusObject import *

class norm(pigasusObject):
    def __init__ ( self, field = None, type = None, func = None, paramevalfunc = False, exact = None ):
        pigasusObject.__init__(self)

        self.id = self.com.nnorms
        self.nparam = 0
        self.paramevalfunc = paramevalfunc

        if field is not None:
            self.field = field
            self.space = field.space
            self.loc_id = self.space.grids.add_norm_id(self)
        else:
            raise("You must give a field for the current norm")

        if type is not None:
            self.type = type
        else:
            self.type = _cst.NORM_L2

        self._set_nparam()

        from .utils import function
        if func is not None:
            self.func = function(func, space=self.space)
        else:
            self.defaultFuncParam()

        if exact is not None:
            self.exact = function(exact, space=self.space)
        else:
            self.defaultFuncExact()

        # this must be the last thing to do
        self.com.nnorms += 1
        self.com.norms.append(self)

    def setInfoData(self):
        """
        prints informations about the current norm
        """
        self.infoData['id'] = str(self.id)
        self.infoData['field'] = str(self.field.id)
        self.infoData['space'] = str(self.space.id)
        self.infoData['loc_id'] = str(self.loc_id)
        self.infoData['nparam'] = str(self.nparam)
        self.infoData['paramevalfunc'] = str(self.paramevalfunc)
        self.infoData['type'] = str(self.type)

    def _getGlobalNorm(self):
        return self.com.pyfem.getglobalnorm ( self.id )

    def _getPatchNorm(self):
        li_npatchs = self.space.grids.npatchs
        return self.com.pyfem._getPatchNorm ( self.id, li_npatchs )

    def _getElementNorm(self, ai_patch):

        li_nel = self.space.grids.list_grid[ai_patch].nel
        return self.com.pyfem._getElementNorm ( self.id, ai_patch, li_nel)

    def get(self, type=0, ai_patch=None):
        """
        returns values for a given type of norm
        type = 0    : for a global computation
        type = 1    : for a patch computation
        type = 2    : for an element computation
        """
        if (type == 0) :
            return self._getGlobalNorm()
        if (type == 1) :
            return self._getPatchNorm()
        if (type == 2) and (ai_patch is not None):
            return self._getElementNorm(ai_patch)

    def setEvalNorm(self, ai_patch=0, fields=[], funcs=[]):
        """
        fields is a list of fields
        funcs is a list of functions
        """
        lpr_pts = self.space.get_points(ai_patch)
        list_pts = []
        for i in range(0, self.space.dim):
            list_pts.append(lpr_pts[i,0,:])
        lpr_pts = list_pts

        li_dim = self.space.dim
        if li_dim not in [2]:
            print("setEvalNorm: Not yet implemetend for the desired dimension")

        lpi_shape = lpr_pts.shape[0:-1]
        lpr_val = _np.zeros((1,lpi_shape[0],lpi_shape[1]))
        for F in fields:
            lpr_f = F.eval(ai_patch, elts)[ai_patch,:,:]
            lpr_val[0,:,:] += lpr_f[:,:]
        for func in funcs:
            lpr_f = _np.zeros(lpr_pts.shape[0:-1])
            for (i,list_p) in enumerate(lpr_pts):
                for (j,p) in enumerate(list_p):
                    lpr_f[i,j] =func (p[0], p[1])[0]
            lpr_val[0,:,:] += lpr_f[:,:]
        self.com.pyfem.set_field_on_grids(self.field.id, ai_patch, lpr_val)

    def _set_nparam(self):

        if ( self.type in [ _cst.NORM_L2 ] ):
            self.nparam = 1
            return
        if ( self.type in [ _cst.NORM_H1 ] ):
            li_dim = self.space.dim
            self.nparam = li_dim**2
            return
        else :
            print("NORM-_set_nparam : type not implemented yet")
            import sys; sys.exit(1)

    def evalfunc(self, ai_patch, apr_points, elts=None, type="param"):
        """
        Evaluation of the param-function over a given list of points
        """
        if not self.paramevalfunc :
            lpr_val = self._evalfunc_std(ai_patch, apr_points, elts, type)
        else:
            lpr_parampts = self.space.get_parametricPoints(ai_patch_id=ai_patch)
            lpr_val = self._evalfunc_std(ai_patch, lpr_parampts, elts, type)
        return lpr_val

    def _evalfunc_std(self, ai_patch, apr_points, elts, type):
        """
        sequential version of the evaluation
        """
        if type == "param":
#            print "==== param evaluation"
            return self.func(apr_points)
        if type == "exact":
#            print "==== exact evaluation"
            return self.exact(apr_points)

    def defaultFuncParam(self):
        li_dim = self.space.dim

        if ( self.type in [ _cst.NORM_L2 ] ):
            if li_dim == 1:
                func = lambda x : [1.0]
            if li_dim == 2:
                func = lambda x,y : [1.0]
            if li_dim == 3:
                func = lambda x,y,z : [1.0]
        elif ( self.type in [ _cst.NORM_H1 ] ):
            if li_dim == 1:
                func = lambda x : [1.0]
            if li_dim == 2:
                func = lambda x,y : [1.0, 0.0, 0.0, 1.0]
            if li_dim == 3:
                func = lambda x,y,z : [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
        else :
            print("NORM-defaultFuncParam : type not implemented yet")
            import sys; sys.exit(1)

        from .utils import function
        self.func = function(func, space=self.space)

    def defaultFuncExact(self):
        li_dim = self.space.dim

        if li_dim == 1:
            func = lambda x : [0.0] * self.field.ndof
        elif li_dim == 2:
            func = lambda x,y : [0.0] * self.field.ndof
        elif li_dim == 3:
            func = lambda x,y,z : [0.0] * self.field.ndof
        else :
            raise("type not implemented yet")

        from .utils import function
        self.exact = function(exact, space=self.space)


    def set_func(self, exact):
        """
        this sets the param-function of the current field
        """
        from .utils import function
        self.exact = function(exact, space=self.space)

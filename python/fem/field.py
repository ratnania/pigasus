# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['field']
__date__ ="$Jan 11, 2012 3:24:13 PM$"

from . import common_obj as _com
from . import constants as _cst
import numpy as _np
from numpy import zeros, asarray, zeros_like
from .pigasusObject import *
from caid.cad_geometry import cad_nurbs

##############################################################################
#
#       Field class
#
##############################################################################
class field(pigasusObject):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        obj.id = None
        return obj

    def __init__ ( self, as_name='F', func = None, pfunc=None, space = None, type = None \
    , operator = None, field = None, parameval = False, paramevalfunc = False \
    , func_arguments = [] ):
        # each field is identified by an id,
        # which will be given when calling iso.addfield
        # as_name : is the name used for visualization diagnostics
        # func : is the exact solution
        # space : is the discrete space
        # type : is the projection type : L2, ...
        # operator : is the tye of the operator, applied to (a) specific field(s)
        # field : the operande for the operator
        # parameval : True if we want to have the eval-basis on the parametric domain : grad,...
        pigasusObject.__init__(self)

        self.id             = self.com.nfields
        self.name           = as_name
        self.ndof           = 1
        self._size          = 0
        self.func_arguments = func_arguments
        self.type           = type
        self.operator       = operator
        self.operande_id    = 0
        self.nderiv         = 1
        self.parameval      = parameval
        self.paramevalfunc  = paramevalfunc
        self.nparam         = 0
        self.func           = None
        self.pfunc          = None

        if space is None:
            print("field: you must specify the space")
            import sys; sys.exit(0)

        self.space = space
        self._size = space.size
        self.loc_id = space.grids.add_field_id(self)

        if self.type is not None :
            if (self.operator is None) or (self.operator == _cst.IDENTITY) :
                self.operator = _cst.IDENTITY

            if (self.type==_cst.FIELD_OPERATOR) and (self.operator in _cst.list_FIELD_OPERATORS) :
                self.nderiv = space.dim
                self.operande_id = field.id
                self.func_arguments = [field]

        else :
            if self.operator is None :
                # default case
                self.type = _cst.PROJECTION_L2
                self.operator = _cst.IDENTITY
                if space.ndof > 1 :
                    self.operator = _cst.IDENTITY_VECT
            if self.operator in _cst.list_FIELD_OPERATORS :
                # this is a Field-Operator
                self.type = _cst.FIELD_OPERATOR
                self.nderiv = space.dim
                self.operande_id = field.id
                self.func_arguments = [field]


        if func is not None:
            self.set_func(func)
        else:
            func = self._set_default_func()
            self.set_func(func)

        if pfunc is not None:
            self.set_pfunc(pfunc)
        else:
            pfunc = self._set_default_func()
            self.set_pfunc(pfunc)

        self._set_nparam()
        self._set_ndof()

        # update the ndof of the corresponding grids
        if ( self.space.grids.ndof < self.ndof ) :
            self.space.grids.ndof = self.ndof
        # this must be the last thing to do
        self.com.nfields += 1
        self.com.fields.append(self)

    @property
    def size(self):
        if self.id is None:
            return len(self._values)
        else:
            return self._size

    def set_func(self, func):
        """
        this sets the param-function of the current field
        """
        li_dim = self.space.dim
        from . import func_tools as ft
        if (self.operator in _cst.list_FIELD_OPERATORS) or (len(self.func_arguments) == 0) :
            from .utils import function
            self.func = function(func, space=self.space)
        else:
            raise("Not used anymore. Dead code")

    def set_pfunc(self, func):
        """
        this sets the param-function of the current field
        """
        if (self.operator in _cst.list_FIELD_OPERATORS) or (len(self.func_arguments) == 0) :
            from .utils import function
            self.pfunc = function(func, space=self.space)
        else:
            raise("Not used anymore. Dead code")

    def setInfoData(self):
        """
        prints informations about the current field
        """
        self.infoData['id'] = str(self.id)
        self.infoData['name'] = str(self.name)
        self.infoData['space'] = str(self.space.id)
        self.infoData['size'] = str(self.size)
        self.infoData['ndof'] = str(self.ndof)
        self.infoData['type'] = str(self.type)
        self.infoData['operator'] = str(self.operator)
        self.infoData['operande_id'] = str(self.operande_id)
        self.infoData['func_arguments'] = str(self.func_arguments)
        self.infoData['nderiv'] = str(self.nderiv)
        self.infoData['parameval'] = str(self.parameval)
        self.infoData['paramevalfunc'] = str(self.paramevalfunc)
        self.infoData['loc_id'] = str(self.loc_id)

    def _set_default_func(self):
        """
        sets the number of parameters for the current field
        """
        dim = self.space.dim
        if dim == 1:
            func = lambda x : [ 1. ]
        if dim == 2:
            func = lambda x,y : [ 1. ]
        if dim == 3:
            func = lambda x,y,z : [ 1. ]
        return func

    def _set_nparam(self):
        self.nparam = self._compute_nparam(self.type, self.operator)

    def _set_ndof(self):
        self.ndof = self._compute_ndof(self.type, self.operator)

    def _compute_nparam(self, ai_type, ai_operator):
        li_nparam = 0
        li_dim = self.space.dim

        if (ai_type == _cst.PROJECTION_L2):
            if ai_operator in [_cst.IDENTITY]:
                li_nparam = 1
            elif ai_operator in [_cst.GRAD, _cst.CURL]:
                li_nparam = li_dim
            else:
                print("field: nparam not yet defined")
                print("for data : ", ai_type, ai_operator)

        if (ai_type == _cst.FIELD_OPERATOR):

            if li_dim  == 1:
                if ai_operator in [_cst.GRAD_S, _cst.SECOND_DERIV_S_FIELD, _cst.SECOND_DERIV_FIELD]:
                    li_nparam = 1
            if li_dim  == 2:
                if ai_operator in [_cst.HESSIAN_FIELD]:
                    li_nparam = 1
                if ai_operator in [_cst.GRAD_S]:
                    li_nparam = 2
                if ai_operator in [_cst.SECOND_DERIV_S_FIELD, _cst.SECOND_DERIV_FIELD]:
                    li_nparam = 3

        return li_nparam

    def _compute_ndof(self, ai_type, ai_operator):
        """
        for scalar fields, dof is equal to 1
        """
        li_dof = 1

        return li_dof

    def get(self):
        """
        returns the coefficients iver the descrite space
        """
        if self.id is None:
            return self._values
        else:
            return self.com.pyfem.getfield ( self.id , self.size )

    def set(self, other):
        """
        sets the current field-coeff to a given value or array
        """
        if self.id is None:
            self._values = other
        else:
            if type(other) is float:
                self._set_val(other)
            if type(other) is _np.ndarray:
                self._set_array(other)
            if _com.isField(other):
                self._set_array(other.get())

    def _set_val(self, ar_val, type="coefs"):
        if type=="coefs":
            self.com.pyfem.field_setval_coefs(self.id, ar_val )
        if type=="values":
            self.com.pyfem.field_setval_values(self.id, ar_val )

    def _set_array(self, apr_val, type="coefs"):
        if type=="coefs":
            self.com.pyfem.field_setarray_coefs(self.id, apr_val )
        if type=="values":
            self.com.pyfem.field_setarray_values(self.id, apr_val )

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
        pvalues = self.pfunc(V.get_sites())
        n,m = values.shape
        v = _np.zeros((1,m))
        for i in range(0,n):
            v += values[i,:] * pvalues[i,:]

        return v

    def varfunc(self, fields=[]):
        """
        this routine defines variables for the field function,
        variables must be of the form  (fields=[...])
        this is to implement functions of the form:
        F(u,v,[x,y]) where u and v are fields
        """
        self.func_arguments = fields

    def _evalfunc_generic(self, ai_patch, apr_points, elts=None):
        # using fields
        if self.type == _cst.FIELD_OPERATOR :
#            print "_evalfunc_generic for OPERATOR"
            if (self.operator in _cst.list_FIELD_OPERATORS):
                return self.eval(ai_patch, elts)
        else :
#            print "Field : _evalfunc_generic for NON LINEAR evaluation"
            lpr_val = self._evalfunc_nonlin(ai_patch, apr_points,elts)
            return lpr_val

    def _evalfunc_nonlin(self, ai_patch, apr_points, elts=None):
        """
        sequential version of the evaluation
        """
        # loop over fields involved in the function
#        li_nfields = len(self.func_arguments)
#        print "li_nfields=", li_nfields

        list_values_F = []
        for F in self.func_arguments:
            list_val = F.eval(ai_patch, elts)
#            print "argfunc id : ", F.id
#            print "argfunc coefs : ", F.get()
#            print "eval on grids ", list_val

            list_values_F.append(list_val)

        # TODO to change when passed to ndof > 1
        lpr_val = self.func(list_values_F, apr_points)
#        print "current Field : ", self.id
#        print "lpr_val=", lpr_val
#        print "lpr_val.shape = ", lpr_val.shape
        return lpr_val

    def set_info(self, patch_id=0, operator=None):
        """
        set new field's info
        this saves new info into fortran
        """
        li_patch = patch_id
        Gr = self.space.grids
        grids_id = Gr.id
        # ...
        # first we change the field's type and operator
        # ...
        if self.type == _cst.FIELD_OPERATOR:
            li_type = self.type
            li_operande_id = self.operande_id
            li_operator = self.operator
            li_nparam = self.nparam
        else:
            li_type = _cst.FIELD_OPERATOR
            li_operande_id = self.id
            if operator is None:
                li_operator = _cst.IDENTITY
            else:
                li_operator = operator
            li_nparam = self._compute_nparam(li_type, li_operator)

        self.type = li_type
        self.operator = li_operator

        li_parameval = 0
        if self.parameval :
            li_parameval = 1
        self.com.pyfem.set_field(self.id, self.ndof, self.space.id, self.size, self.loc_id \
        , li_type, li_operator, li_operande_id, li_parameval, li_nparam)
        # ...

    def evaluate(self, patch_id=None, sites=None, elts=None, fmt=True, nderiv=0,
                 operator=None, parametric=False):
#        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        V   = self.space
        geo = V.geometry
#        print geo
        G   = V.grids

        if operator == "grad":
            nderiv=1
        if operator == "second_deriv":
            nderiv=2

        if patch_id is None:
            patch = V.currentPatch
            patch_id = V.geometry.index(patch)
#        print "patch_id : ", patch_id

        if sites is not None:
            _sites = sites
        else:
            _sites = G.get_sites()[patch_id]

        npts = G.k

        nrb = geo[patch_id]
        nx = G.list_grid[patch_id].nx

        knots   = nrb.knots
        weights = nrb.weights
        shape = list(nrb.shape)
        _C = self.tomatrix(patch_id)
        C = zeros(shape+[3])
        C[...,0] = _C
        srf = cad_nurbs(knots, C, weights=weights)

        # ... defining the evaluate function
        APPEND = False
        NOUT   = 1
        if nderiv == 0:
            EVALUATE = srf.__call__
            APPEND = True
            NOUT   = 1
        if nderiv > 0:
            def EVALUATE (u=None, v=None, w=None):
                return srf.evaluate_deriv(u=u, v=v, w=w, nderiv=nderiv)
        # ...
#        print "======================="
#        print "INPUT"

        from .utils import evaluator
        D = evaluator(nrb.dim, _sites, EVALUATE, APPEND, nx=nx, npts=npts,
                         elts=elts, fmt=fmt)
#        print "OUTPUT"
#        print D[0,0,:]
#        print D.shape
#        print "======================="

        ndim = len(shape)

        xyz_arrays = None
        if operator is None:
            if APPEND:
                xyz_arrays = D[0,0]
            else:
                xyz_arrays = D
        else:

            if operator == "grad":
                if ndim == 1:
                    nb = 1
                    ne = 2
                if ndim == 2:
                    nb = 1
                    ne = 3
                if ndim == 3:
                    nb = 1
                    ne = 4

            if operator == "second_deriv":
                if ndim == 1:
                    nb = 2
                    ne = 3
                if ndim == 2:
                    nb = 3
                    ne = 6
                if ndim == 3:
                    nb = 4
                    ne = 10

            xyz_arrays = []
            for i in range(nb,ne):
                du = D[0,i,:]
                xyz_arrays.append(du)

        # returns values on the logical domain
        if parametric:
            return xyz_arrays
        # returns values on the physical domain
        else:
            # ...
            P = V.get_points()

            # ...
            if nderiv >= 1:
                xyz    = P[:,0,:]
                U      = xyz_arrays[0,0,:]

                if ndim == 1:
                    xyzdu  = P[:,1,:]
                    Udu    = xyz_arrays[0,1,:]
                    print("Physical evaluation not yet implemented")

                if ndim == 2:
                    xyzdu  = P[:,1,:]
                    xyzdv  = P[:,2,:]

                    Udu    = xyz_arrays[0,1,:]
                    Udv    = xyz_arrays[0,2,:]

                    xdu    = xyzdu[0,...]
                    ydu    = xyzdu[1,...]
                    xdv    = xyzdv[0,...]
                    ydv    = xyzdv[1,...]

                    jac = xdu * ydv - xdv * ydu

                    Udx =   ydv * Udu - ydu * Udv
                    Udx /= jac
                    Udy = - xdv * Udu + xdu * Udv
                    Udy /= jac

                if ndim == 3:
                    xyzdu  = P[:,1,:]
                    xyzdv  = P[:,2,:]
                    xyzdw  = P[:,3,:]

                    Udu    = xyz_arrays[0,1,:]
                    Udv    = xyz_arrays[0,2,:]
                    Udw    = xyz_arrays[0,3,:]

                    print("Physical evaluation not yet implemented")
            # ...

            # ...
            if nderiv >= 2:
                if ndim == 1:
                    xyzduu  = P[:,2,:]

                    Uduu    = xyz_arrays[0,2,:]
                if ndim == 2:
                    xyzduu  = P[:,3,:]
                    xyzduv  = P[:,4,:]
                    xyzdvv  = P[:,5,:]

                    Uduu    = xyz_arrays[0,3,:]
                    Uduv    = xyz_arrays[0,4,:]
                    Udvv    = xyz_arrays[0,5,:]

                    xdu    = xyzdu[0,...]
                    ydu    = xyzdu[1,...]
                    xdv    = xyzdv[0,...]
                    ydv    = xyzdv[1,...]
                    xduu   = xyzduu[0,...]
                    yduu   = xyzduu[1,...]
                    xduv   = xyzduv[0,...]
                    yduv   = xyzduv[1,...]
                    xdvv   = xyzdvv[0,...]
                    ydvv   = xyzdvv[1,...]

                    C1 = Uduu - xduu * Udx - yduu * Udy
                    C2 = Uduv - xduv * Udx - yduv * Udy
                    C3 = Udvv - xdvv * Udx - ydvv * Udy
                    Udxx =   C1 * ydv**2    - 2 * C2 * ydu * ydv + C3 * ydu**2
                    Udxx /= jac**2
                    Udxy = - C1 * xdv * ydv + C2 *(xdu * ydv + xdv * ydu) - C3 * xdu * ydu
                    Udxy /= jac**2
                    Udyy =   C1 * xdv**2    - 2 * C2 * xdu * xdv + C3 * xdu**2
                    Udyy /= jac**2

                if ndim == 3:
                    xyzduu  = P[:,4,:]
                    xyzdvv  = P[:,5,:]
                    xyzdww  = P[:,6,:]
                    xyzduv  = P[:,7,:]
                    xyzdvw  = P[:,8,:]
                    xyzduw  = P[:,9,:]

                    Uduu    = xyz_arrays[0,4,:]
                    Udvv    = xyz_arrays[0,5,:]
                    Udww    = xyz_arrays[0,6,:]
                    Uduv    = xyz_arrays[0,7,:]
                    Udvw    = xyz_arrays[0,8,:]
                    Uduw    = xyz_arrays[0,9,:]
            # ...

            # ... copy data into xyz_arrays
            if nderiv >= 1:
                xyz_arrays[0,0,:] = U

                if ndim == 1:
                    xyz_arrays[0,1,:] = Udx

                if ndim == 2:
                    xyz_arrays[0,1,:] = Udx
                    xyz_arrays[0,2,:] = Udy

                if ndim == 3:
                    xyz_arrays[0,1,:] = Udx
                    xyz_arrays[0,2,:] = Udy
                    xyz_arrays[0,3,:] = Udz
            # ...
            # ...
            if nderiv >= 2:
                if ndim == 1:
                    xyz_arrays[0,2,:] = Udxx

                if ndim == 2:
                    xyz_arrays[0,3,:] = Udxx
                    xyz_arrays[0,4,:] = Udxy
                    xyz_arrays[0,5,:] = Udyy

                if ndim == 3:
                    xyz_arrays[0,4,:] = Udxx
                    xyz_arrays[0,5,:] = Udyy
                    xyz_arrays[0,6,:] = Udzz
                    xyz_arrays[0,7,:] = Udxy
                    xyz_arrays[0,8,:] = Udyz
                    xyz_arrays[0,9,:] = Udxz
            # ...

            return xyz_arrays

    def grad(self, patch_id=None, apr_points=None, elts=None):
        gradF = self.evaluate(patch_id=patch_id, elts=elts, operator="grad")
        return gradF

    def fast_plot(self, ai_patch_id = 0, useControlPoints=True,
                  savedPoints=True, list_P=None, withpcolor=False, vmin=None, vmax=None):
        """
        use this function for a fast plot of the current field.
        This affects the coeffs to either the control points, or the averages knots images.
        """

        if self.space.dim==1:
            self._fast_plot_1D(ai_patch_id = ai_patch_id \
            , useControlPoints=useControlPoints\
            , savedPoints=savedPoints\
            , list_P=list_P)

        if self.space.dim==2:
            self._fast_plot_2D(ai_patch_id = ai_patch_id \
                               , useControlPoints=useControlPoints\
                               , savedPoints=savedPoints\
                               , list_P=list_P \
                              , withpcolor=withpcolor \
                              , vmin=vmin \
                              , vmax=vmax)

    def _fast_plot_1D(self, ai_patch_id = 0, useControlPoints=True, savedPoints=True, list_P=None):
        """
        input :
            ai_patch_id : the patch id
            useControlPoints : use the control points as an approximation of the evaluation points
            savedPoints : used if useControlPoints is False.
                savedPoints = True => uselist_P as an input for the evaluation grid
                savedPoints = False => must compute the evaluation grid


        X is the 1D grid
        """
        # ...
        def knot_average(nrb):
            list_Gr = []
            for li_d in range(0,nrb.dim):
                list_knots = nrb.knots[li_d]
                li_N = nrb.shape[li_d]
                li_P = nrb.degree[li_d]

                list_t = []
                for i in range(1, li_N+1):
                    lr_t = sum(list_knots[i:i+li_P]) / li_P
                    list_t.append(lr_t)
                list_Gr.append(list_t)
            return list_Gr
        # ...

        from pylab import plot
        V = self.space
        # ...
        # getting the average knots and their images
        # ...
        li_patch_id = ai_patch_id
        nrb = V.geometry[li_patch_id]
        list_xi = knot_average(nrb)

        lpi_n = nrb.shape

        # ...
        # using Control Points for the correspondance
        # ...
        if useControlPoints:
            lpi_shape = list(lpi_n)+[nrb.points.shape[-1]]
            lpr_P = _np.asarray(nrb.points).reshape(lpi_shape).transpose()
            X = lpr_P[0,:]
        else :
            print("Warning fast_plot: not yet implemented with the exact evaluation points, ie average knots images")
        # ...

        list_faces_dir = V.list_faces_dir
        list_faces_duplicated = V.list_faces_duplicated
        list_faces_duplicata  = V.list_faces_duplicata

        if len(list_faces_duplicated) > 0 :
            print("Warning fast_plot: not yet implemented with duplicated faces boundary conditions")

        lpr_val = _np.zeros(lpi_n)

        lpi_m = _np.asarray(lpi_n)
        list_f_dir = list_faces_dir[li_patch_id]
        list_f_duplicated = []
        list_f_duplicata  = []

        nx_min =  0 ; nx_max = lpi_n[0]

        if 0 in list_f_dir:
            lpi_m[0] -= 1
            nx_min += 1
        if 1 in list_f_dir:
            lpi_m[0] -= 1
            nx_max -= 1

        lpr_u = self.get().reshape(lpi_m)

        for i in range (nx_min, nx_max):
            i_ = i
            if nx_min == 1:
                i_ = i - 1
            lpr_val [i] = lpr_u[i_]

        # ...
        #G = _np.meshgrid(_np.asarray(list_xi), _np.asarray(list_eta) )
        #X = G[0] ; Y = G[1]
        #pl.plot(X,Y,'o')
        plot(X, lpr_val)
        # ...

    def _fast_plot_2D(self, ai_patch_id = 0, useControlPoints=True,
                      savedPoints=True, list_P=None, withpcolor=False,
                      vmin=None, vmax=None):
        """
        input :
            ai_patch_id : the patch id
            useControlPoints : use the control points as an approximation of the evaluation points
            savedPoints : used if useControlPoints is False.
                savedPoints = True => uselist_P as an input for the evaluation grid
                savedPoints = False => must compute the evaluation grid


        X,Y is the 2D grid
        """
        # ...
        def knot_average(nrb):
            list_Gr = []
            for li_d in range(0,nrb.dim):
                list_knots = nrb.knots[li_d]
                li_N = nrb.shape[li_d]
                li_P = nrb.degree[li_d]

                list_t = []
                for i in range(1, li_N+1):
                    lr_t = sum(list_knots[i:i+li_P]) / li_P
                    list_t.append(lr_t)
                list_Gr.append(list_t)
            return list_Gr
        # ...

        V = self.space
        # ...
        # getting the average knots and their images
        # ...
        li_patch_id = ai_patch_id
        nrb = V.geometry[li_patch_id]
        list_xi, list_eta = knot_average(nrb)

        lpi_n = nrb.shape

        # ...
        # using Control Points for the correspondance
        # ...
        if useControlPoints:
            X = nrb.points[:,:,0]
            Y = nrb.points[:,:,1]
        else :
            print("Warning fast_plot: not yet implemented with the exact evaluation points, ie average knots images")
        # ...

        lpr_val = self.tomatrix(li_patch_id)
#        lpr_val = V.inv_vect(0, self.get())
#        _np.savetxt('u.txt', lpr_val)
#        _np.savetxt('X.txt', X)
#        _np.savetxt('Y.txt', Y)

        # ...
        from pylab import contourf, pcolor
        if not withpcolor:
            contourf(X, Y, lpr_val)
        else:
            pcolor(X, Y, lpr_val, vmin=vmin, vmax=vmax)
        # ...

    def plot(self, withpcolor=False, n=None):
        from matplotlib.pyplot import contourf, pcolor
        from numpy import min, max, linspace, zeros_like
        geo = self.space.geometry
        U = self
        vmin = min(U.get())
        vmax = max(U.get())
        dim = geo.dim
        for patch_id in range(0, geo.npatchs):
            nrb = geo[patch_id]
            list_n = n
            if n is None:
                list_n = [50]*nrb.dim

            C = zeros_like(nrb.points)

            x = nrb.points[...,0]
            y = nrb.points[...,1]
            z = U.tomatrix(patch_id) #.transpose()
            C[...,0] = x
            C[...,1] = y
            C[...,2] = z
            srf = cad_nurbs(nrb.knots, C, weights=nrb.weights)
            list_t = []
            for axis in range(0, nrb.dim):
                t = linspace(nrb.knots[axis][0],nrb.knots[axis][-1],list_n[axis])
                list_t.append(t)
            P = srf(u=list_t[0], v=list_t[1])

            if dim == 1:
                x = P[...,0]
                y = P[...,1]
                plot(x,y)

            if dim == 2:
                x = P[...,0]
                y = P[...,1]
                z = P[...,2]
                if withpcolor:
                    pcolor(x,y,z,vmin=vmin,vmax=vmax)
                else:
                    contourf(x,y,z)

    def reset(self, coef=True, values=False):
        """
        this sets the field coeffs or values to 0
        """

        if coef :
            self.com.pyfem.pyfem_reset_field(self.id)
        if values :
            print("Field reset : not implemente yet")

    def save(self, filename):
        _np.savetxt(filename, self.get())

    def load(self, filename):
        self.set(_np.genfromtxt(filename))

    def export(self, filename):
        from caid.cad_geometry import cad_geometry, cad_nurbs
        from caid.field import field as iga_field
        U = self
        geo   = U.space.geometry
        geo_c = cad_geometry()
        for i in range(0, geo.npatchs):
            nrb = geo[i]
            C = zeros_like(nrb.points)
            _C = U.tomatrix(i)
            shape = list(nrb.shape)
            C = zeros(shape+[3])
            C[...,0] = _C
            srf = cad_nurbs(nrb.knots, C, weights= nrb.weights)
            geo_c.append(srf)
        F = iga_field(geo, geo_c)
        F.save(filename)

    def tomatrix(self, ai_patch):
        """
        transforms the coeff 1D-array to a 2D-matrix
        """
        V = self.space
#        print "------------"
#        print "geo.npatchs : ", V.geometry.npatchs
#        print "patch id : ", ai_patch
#        print "dim : ", V.dim
#        print "shape : ", V.geometry[ai_patch].shape
        if V.dim == 1 :
            [li_n_1] = V.geometry[ai_patch].shape
            return self.com.pyfem.field_to_matrix_1d ( self.id, ai_patch  \
            , li_n_1 )
        if V.dim == 2 :
            [li_n_1, li_n_2] = V.geometry[ai_patch].shape
            return self.com.pyfem.field_to_matrix_2d ( self.id, ai_patch  \
            , li_n_1, li_n_2 )
        if V.dim == 3 :
            [li_n_1, li_n_2, li_n_3] = V.geometry[ai_patch].shape
            return self.com.pyfem.field_to_matrix_3d ( self.id  \
            , ai_patch, li_n_1, li_n_2, li_n_3 )

    def frommatrix(self, ai_patch, C, alpha=1.0):
        """
        transforms the coeff nD-matrix to 1D-array
        """
        V = self.space
        if V.dim == 1 :
            [li_n_1] = V.geometry[ai_patch].shape
            self.com.pyfem.field_from_matrix_1d ( self.id, ai_patch \
                                                 , alpha \
                                                 , C   \
                                                 , li_n_1 )
        if V.dim == 2 :
            [li_n_1, li_n_2] = V.geometry[ai_patch].shape
            self.com.pyfem.field_from_matrix_2d ( self.id, ai_patch \
                                                 , alpha \
                                                 , C  \
                                                 , li_n_1, li_n_2 )
        if V.dim == 3 :
            [li_n_1, li_n_2, li_n_3] = V.geometry[ai_patch].shape
            self.com.pyfem.field_from_matrix_3d ( self.id, ai_patch \
                                                 , alpha \
                                                 , C \
                                                 , li_n_1, li_n_2, li_n_3 )

    def set_initialize(self):
        self.com.pyfem.pyfem_set_initialize_field(self.id)

    def set_finalize(self):
        self.com.pyfem.pyfem_set_finalize_field(self.id)

    def tocsr(self):
        from scipy import sparse

        n = self.size
        m = 1
        nnz = n

        I = array(list(range(0,n)))
        J = array([0]*nnz)
        V = self.get()

        return sparse.coo_matrix((V,(I,J)),shape=(n,m)).tocsr()

    def riseFromBoundary(self, patch_id, list_faces, g=None, list_g=None):
        from . import bnd_function as bf
        V = self.space
        nrb = V.geometry[patch_id]
        size = asarray(nrb.shape).prod()

        patch_id = self.space.geometry.index(nrb)
#        print "************"
#        print "patch ", patch_id
#        print "faces ", list_faces
#        saved = self.get()
#        self.reset()

#        lpr_g = _np.zeros(V.size)

        nfaces = len(list_faces)
        if g is not None:
            _list_g = [g] * nfaces
        elif list_g is not None:
            _list_g = list_g
        assert(len(_list_g)==len(list_faces))

        # ...
        def _update_coeff(face, done_faces, gh):
            # nothing to do if dim==1

            # must transpose before treatment
            gh = gh.transpose()

            if V.dim == 2 :
                if (face == 0):
                    if done_faces[1]:
                        gh[0,0]   = 0.
                    if done_faces[3]:
                        gh[-1,0]  = 0.

                if (face == 1):
                    if done_faces[0]:
                        gh[0,0]   = 0.
                    if done_faces[2]:
                        gh[0,-1]  = 0.

                if (face == 2):
                    if done_faces[1]:
                        gh[0,-1]  = 0.
                    if done_faces[3]:
                        gh[-1,-1] = 0.

                if (face == 3):
                    if done_faces[0]:
                        gh[-1,0]  = 0.
                    if done_faces[2]:
                        gh[-1,-1] = 0.

            # re-transpose to get the right array
            gh = gh.transpose()
            return gh
        # ...

        # ...
        if V.dim == 1:
            done_faces = [False] * 2
        if V.dim == 2:
            done_faces = [False] * 4
        if V.dim == 3:
            done_faces = [False] * 6

        mat_g = zeros(nrb.shape)
        for (face, f) in zip(list_faces, list_g):
            bnd_fct = bf.bnd_function(nrb, face, f)
            gh      = bnd_fct.rise("coeff")
            done_faces[face] = True
            gh = _update_coeff(face, done_faces,gh)
            gh = gh.transpose()
            mat_g += gh
#            self.reset()
#            self.frommatrix(patch_id,gh)
#            lpr_g += self.get()
        # ...
        return mat_g

    def __call__(self, u=None, v=None, w=None, patch_id=None):
        geo = self.space.geometry
        U = self
        dim = geo.dim
        uvw = [u,v,w][:dim]

        list_patchs = []
        if patch_id is None:
            list_patchs = [nrb for nrb in geo]
        else:
            list_patchs = [geo[patch_id]]

        list_P = []
        for nrb in list_patchs:
            nrb_id = geo.index(nrb)

            C = zeros_like(nrb.points)
            x = nrb.points[...,0]
            y = nrb.points[...,1]
            z = U.tomatrix(nrb_id) #.transpose()
            C[...,0] = x
            C[...,1] = y
            C[...,2] = z
            srf = cad_nurbs(nrb.knots, C, weights=nrb.weights)
            P = srf(u=u,v=v,w=w)
            if patch_id is not None:
                return P[...,-1]
            else:
                list_P.append(P[...,-1])
        return list_P

    def __pos__(self):
        return self.__mul__(1.0)

    def __neg__(self):
        return self.__mul__(-1.0)

    def __sub__(self, other):
        return self.__add__(-other)

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __isub__(self, other):
        return self.__iadd__(-other)

    def __add__(self, other):
        F = field.__new__(type(self))
        if _com.isFloat(other):
            F.set(self.com.pyfem.field_scaladdition2(self.id, other, self.size ))
        elif (self.id is None) or (other.id is None):
            F.set(self.get() + other.get())
        elif _com.isField(other):
            F.set(self.com.pyfem.field_addition2(self.id, other.id, self.size ))
        return F

    def __mul__(self, other):
        F = field.__new__(type(self))
        if _com.isFloat(other):
            F.set(other * self.get())
        elif (self.id is None) or (other.id is None):
            F.set(self.get() * other.get())
        elif _com.isField(other):
            F.set(self.com.pyfem.field_multiplication2(self.id, other.id, self.size))
        return F

    def __iadd__(self, other):
        if _com.isFloat(other):
            self.com.pyfem.field_add_scal(self.id, other, self.id )
        if _com.isField(other):
            self.com.pyfem.field_addition(self.id, other.id, self.id )
        return self

    def __imul__(self, other):
        if _com.isFloat(other):
            self.com.pyfem.field_mult_scal(self.id, other, self.id)
        if _com.isNumpyArray(other):
            self.com.pyfem.field_mult_array (self.id, other)
        if _com.isField(other):
            self.com.pyfem.field_mult_field(self.id, other.id, self.id )
        if _com.isMatrix(other):
            self.com.pyfem.matrix_mult_field(other.id,self.id,self.id)
        return self
##############################################################################


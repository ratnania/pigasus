# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['metric']
__date__ ="$Jul 24, 2012 4:22:12 PM$"

from . import common_obj as _com
import numpy as np
from .pigasusObject import *
import sys
from numpy import zeros, asarray

class metric(pigasusObject):
    def __init__(self, analytic=None, points=None, geometry=None, with_igakit=True):
        """
        a metric can be constructed either by giving
        - the analytic form [F,J_F]
        where F is the geometric function and J_F its jacobian
        - an array of points of the form Points[i, 0:Rd,0:Rd], i is the i^th point
            Points[i,0,:] is the physical point
            Points[i,1,:] is the x-derivative
            Points[i,2,:] is the y-derivative
            Points[i,3,:] is the z-derivative
        - a cad_geometry object
        """
        pigasusObject.__init__(self)

        self.type = None
        self.geometry = geometry
        self.F = None
        self.DF = None

        if analytic is not None :
            self.type = "analytic"
#            li_dim = 2
#            import func_tools as ft
#            self.F = ft.internal_func (analytic[0], li_dim)
#            self.DF = ft.internal_func (analytic[1], li_dim)
            self.F = analytic[0]
            self.DF = analytic[1]

        if points is not None :
            self.type = "points"
            self.points = points

        if geometry is not None:
            if with_igakit:
                self.type = "cad_geometry"
            else:
                self.type = "cad_mapping"
                from . import mapping as mp
                self.mapping = mp.mapping(geometry)

        # this must be the last thing to do
        self.id = self.com.nmetrics
        self.com.nmetrics += 1
        self.com.metrics.append(self)

#    def update(self, ai_patch_id, sites, nel, npts):
#        """
#        updates the points and the jacobian for a given patch/sub domain
#        """
#        # apr_points do not have the same shape depending on self.type
#        if self.type == "analytic" :
#            return self._update_analytic(ai_patch_id,sites, nel, npts)
#        if self.type == "points" :
#            return self._update_points(ai_patch_id,sites, nel, npts)

    def update_analytic(self, ai_patch_id, sites, nel, npts):
        xyz = self.F(*sites)
        D = self.DF(*sites)

#        print "--xyz--"
#        print xyz
#        print len(xyz), xyz[0].shape
#        print "---D---"
#        print D
#        print len(D), D[0].shape


        lpr_F = list(zip(*xyz))
        lpr_DF = list(zip(*D))

        lpr_F = np.asarray([np.asarray(list(F)) for F in lpr_F])
        lpr_DF = np.asarray([np.asarray(list(DF)) for DF in lpr_DF])

#        print "lpr_F.shape = ", lpr_F.shape
#        print "lpr_DF.shape = ", lpr_DF.shape

        lpi_shape = [nel, npts]

        list_shape = list(lpi_shape)

        # ...
        # lpr_F treatment
        # ...
        li_n = lpr_F.shape[1]
        lpr_val = np.zeros(list_shape+[li_n])
        for i in range(0,li_n):
                lpr_tmp = lpr_F[:,i]
                lpr_tmp = lpr_tmp.reshape(list_shape)
                lpr_val[:,:,i] = lpr_tmp
        lpr_F = lpr_val
        # ...

        # ...
        # lpr_DF treatment
        # ...
        li_n = lpr_DF.shape[1]
        lpr_val = np.zeros(list_shape+[li_n])
        for i in range(0,li_n):
                lpr_tmp = lpr_DF[:,i]
                lpr_tmp = lpr_tmp.reshape(list_shape)
                lpr_val[:,:,i] = lpr_tmp
        lpr_DF = lpr_val
        # ...

        self.com.pyfem.set_metric_points(self.id, ai_patch_id, lpr_F, lpr_DF)

        x = xyz[0] ; y = xyz[1]
        xdu = D[0] ; xdv = D[1] ; ydu = D[2] ; ydv = D[3]
        lpr_pts = zeros((3,3,len(x)))
        lpr_pts[0,0,:] = x
        lpr_pts[1,0,:] = y
        lpr_pts[0,1,:] = xdu
        lpr_pts[1,1,:] = ydu
        lpr_pts[0,2,:] = xdv
        lpr_pts[1,2,:] = ydv

        return lpr_pts

    def update_cad_geometry(self, ai_patch_id, sites, nderiv, k, nx, nel, npts):
        srf = self.geometry[ai_patch_id]

        # ... defining the evaluate function
#        EVALUATE = srf.evaluate_deriv
        def EVALUATE(u=None, v=None, w=None):
            return srf.evaluate_deriv(u=u, v=v, w=w \
                       , fields=None, nderiv=nderiv, rationalize=srf.rational)
#        APPEND = True
        APPEND = False
        # ...

        from .utils import evaluator
        lpr_pts = evaluator(srf.dim, sites, EVALUATE, APPEND, nx=nx, npts=k,
                      fmt=True)

        x   = lpr_pts[0,0,:]
        y   = lpr_pts[1,0,:]
        xdu = lpr_pts[0,1,:]
        ydu = lpr_pts[1,1,:]
        xdv = lpr_pts[0,2,:]
        ydv = lpr_pts[1,2,:]

        xyz = [x,y]
        D = [xdu,xdv,ydu,ydv]

        lpr_F = list(zip(*xyz))
        lpr_DF = list(zip(*D))

        lpr_F = asarray([asarray(list(F)) for F in lpr_F])
        lpr_DF = asarray([asarray(list(DF)) for DF in lpr_DF])

        lpi_shape = [nel, npts]

        list_shape = list(lpi_shape)

        # ...
        # lpr_F treatment
        # ...
        li_n = lpr_F.shape[1]
        lpr_val = zeros(list_shape+[li_n])
        for i in range(0,li_n):
                lpr_tmp = lpr_F[:,i]
                lpr_tmp = lpr_tmp.reshape(list_shape)
                lpr_val[:,:,i] = lpr_tmp
        lpr_F = lpr_val
        # ...

        # ...
        # lpr_DF treatment
        # ...
        li_n = lpr_DF.shape[1]
        lpr_val = zeros(list_shape+[li_n])
        for i in range(0,li_n):
                lpr_tmp = lpr_DF[:,i]
                lpr_tmp = lpr_tmp.reshape(list_shape)
                lpr_val[:,:,i] = lpr_tmp
        lpr_DF = lpr_val
        # ...
        self.com.pyfem.set_metric_points(self.id, ai_patch_id, lpr_F, lpr_DF)

        return lpr_pts


#    def _update_points(self, ai_patch_id=None, apr_points=None):
#        patch_id = 0
#        if ai_patch_id is not None:
#            patch_id = ai_patch_id
#        points = self.points
#        if apr_points is not None:
#            points = apr_points
#        self.com.pyfem.set_metric_points_advanced(self.id, patch_id, points)

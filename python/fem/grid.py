# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['grid']
__date__ ="$Jan 13, 2012 11:48:35 AM$"

from . import common_obj as _com
import numpy as _np
from .pigasusObject import *
from .quadratures import *

class grid(pigasusObject):
    def __init__(self,ai_current_grids, ai_id, ao_patch \
    , api_k=None, tensorlevel=1, as_type="lobatto" \
    , list_nodes=None, list_w=None):
        pigasusObject.__init__(self)

        if api_k is None:
            print("grid : api_k must be given")
            import sys
            sys.exit(2)

        self.current_grids = ai_current_grids

        self.id = ai_id
        self.patch = ao_patch
        self.tensorlevel = tensorlevel
        self.dim = self.patch.dim
        self.Rd = self.patch.points.shape[-1]
#        li_nel, lpi_nel = ao_patch.get_nel()
#        self.nel = li_nel

        self.list_u = []

        self.maxnpts = 0
        self.dirmaxnpts = 0
        self.nx = []
        self.maxnptsx = []
        self.x  = []
        self.wx = []
        self.grid1D = []

#        print "self.tensorlevel = ", self.tensorlevel

        if self.tensorlevel in [1,2]:
            self._gen_nx()
            if list_nodes is None:
                self._gen_maxnpts(api_k=api_k)
                self._gen_grid(api_k, as_type)
            else :
                self._gen_maxnpts(list_nodes=list_nodes)
                self._gen_grid_from_list(list_nodes, list_w)
        else:
            print("grid : Not yet implemented")
            import sys
            sys.exit(2)

        self.nel = 1
        for nx in self.nx:
            self.nel *= nx
        self.npts = 1
        for k in api_k:
            self.npts *= k+1

    def _gen_nx(self):
        if self.tensorlevel == 0:
            print("tensorlevel = 0 not yet imlpemented!!")
            import sys; sys.exit(0)
        else:
            li_d = 0
            li_n = self.patch.shape[li_d]
            li_p = self.patch.degree[li_d]
            li_un = len(_np.unique(self.patch.knots[li_d]))
            ll_isperiodic = (li_un == li_n + li_p + 1)
            if ll_isperiodic :
                list_u = _np.array([self.patch.knots[li_d][i] for i in range(li_p,li_n+1) ])
                li_nel = li_n - li_p
            else:
                list_u = _np.unique(self.patch.knots[li_d])
                li_nel = len(list_u)-1
            self.list_u.append(list_u)

            if self.tensorlevel == 1:
                self.nx.append(li_nel)

            if self.tensorlevel == 2:
                self.nx.append(1)

            for li_d in range (1,self.dim):
                li_n = self.patch.shape[li_d]
                li_p = self.patch.degree[li_d]
                li_un = len(_np.unique(self.patch.knots[li_d]))
                ll_isperiodic = (li_un == li_n + li_p + 1)
                if ll_isperiodic :
                    list_u = _np.array([self.patch.knots[li_d][i] for i in range(li_p,li_n+1) ])
                    li_nel = li_n - li_p
                else:
                    list_u = _np.unique(self.patch.knots[li_d])
                    li_nel = len(list_u)-1
                self.nx.append(li_nel)
                self.list_u.append(list_u)

    def _gen_maxnpts(self, api_k=None, list_nodes=None):

        if api_k is not None:
            self.maxnpts = 1

            if self.tensorlevel == 0:
                print("tensorlevel = 0 not yet imlpemented!!")
                import sys; sys.exit(0)
            else:
                if self.tensorlevel == 1:
                    li_d = 0
                    li_maxnptsx = api_k[li_d] + 1
                    self.maxnptsx.append(li_maxnptsx)
                    self.maxnpts *= li_maxnptsx
                    self.dirmaxnpts = li_maxnptsx

                if self.tensorlevel == 2:
                    li_d = 0
                    # ...
                    li_n = self.patch.shape[li_d]
                    li_p = self.patch.degree[li_d]
                    li_un = len(_np.unique(self.patch.knots[li_d]))
                    ll_isperiodic = (li_un == li_n + li_p + 1)
                    if ll_isperiodic :
                        list_u = _np.array([self.patch.knots[li_d][i] for i in range(li_p,li_n+1) ])
                        li_nel = li_n - li_p
                    else:
                        list_u = _np.unique(self.patch.knots[li_d])
                        li_nel = len(list_u)-1
                    # ...
                    li_maxnptsx = (api_k[li_d] + 1) * li_nel
                    self.maxnptsx.append(li_maxnptsx)
                    self.maxnpts *= li_maxnptsx
                    self.dirmaxnpts = li_maxnptsx

                for li_d in range (1,self.dim):
                    li_maxnptsx = api_k[li_d] + 1
                    self.maxnptsx.append(li_maxnptsx)

                    self.maxnpts *= li_maxnptsx
                    self.dirmaxnpts = max(self.dirmaxnpts , li_maxnptsx)

        if list_nodes is not None:
            self.maxnptsx = []
            for li_d in range(0, self.dim):
                self.maxnptsx.append(max([len(L) for L in list_nodes[li_d]]))
            self.dirmaxnpts = max(self.maxnptsx)
            self.maxnpts = _np.asarray(self.maxnptsx).prod()

#            print "self.maxnptsx = ", self.maxnptsx
#            print "self.maxnpts = ", self.maxnpts
#            print "self.dirmaxnpts = ", self.dirmaxnpts


    def _gen_vector(self, ai_d, ai_k, as_type):
        qd = quadratures()
        li_nel = self.nx[ai_d]
        li_maxnpts = self.maxnptsx[ai_d]
        list_u = self.list_u[ai_d]

        lpr_local = _np.zeros((li_nel,li_maxnpts+1), dtype=_np.double)
        lpr_localw = _np.ones((li_nel,li_maxnpts), dtype=_np.double)

#        print "_gen_vector : ai_d, li_nel =", ai_d, li_nel

        [lpr_x,lpr_w] = qd.generate(_np.asarray(list_u), ai_k, as_type)

        if (self.tensorlevel==1) or (self.tensorlevel==2 and ai_d>0):
            for li_i in range(0,li_nel):
                lpr_local [li_i,0] = ai_k + 1
                lpr_local [li_i,1:li_maxnpts+1] = lpr_x[li_i,:]
                lpr_localw[li_i,:] = lpr_w[li_i,:]
        if (self.tensorlevel==2 and ai_d==0):
            # we have one element
            lpr_local [0,0] = (ai_k + 1) * (len(list_u)-1)
            for j in range(0, len(list_u)-1):
#                print j
#                print lpr_x[j,:]
#                print lpr_w[j,:]
                lpr_local [0,1+j*(ai_k+1) : 1+(j+1)*(ai_k+1)] = lpr_x[j,:]
                lpr_localw[0,j*(ai_k+1) : (j+1)*(ai_k+1)] = lpr_w[j,:]

        # supprimer la redondance, en mettant le poid a 0
        if as_type in ["uniform"]:
            for li_i in range(0,li_nel-1):
                lpr_localw[li_i, -1] = 0.0

        return lpr_local, lpr_localw, lpr_x

    def _gen_grid(self, api_k, as_type):
        for li_d in range (0,self.dim):
            li_k = api_k[li_d]
            lpr_localx, lpr_localwx, lpr_x = self._gen_vector(li_d, li_k, as_type)

            self.x.append(lpr_localx)
            self.wx.append(lpr_localwx)
            self.grid1D.append(lpr_x)

    def _gen_grid_from_list(self, list_nodes, list_w):
        for (Lx,Lw) in zip(list_nodes, list_w):
            L_x = [] ; L_w = []
            for (lx,lw) in zip(Lx, Lw):
                L_x.append([len(lx)]+list(lx))
                L_w.append(lw)
            self.x.append(L_x)
            self.wx.append(L_w)

    def save_grid(self):
        if self.dim==1:
            self._save_grid_1D()
        if self.dim==2:
            self._save_grid_2D()
        if self.dim==3:
            self._save_grid_3D()
        if self.dim not in [1,2,3]:
            print("save_grid : dimension not done yet")

    def _save_grid_1D(self):
        self.com.pyfem.pyfem_create_from_tensor_product_1d(self.current_grids \
        , self.id \
        , self.x[0], self.wx[0] \
        , self.nx[0], self.maxnptsx[0]+1)

    def _save_grid_2D(self):
        self.com.pyfem.pyfem_create_from_tensor_product_2d(self.current_grids \
        , self.id \
        , self.x[0], self.wx[0] \
        , self.x[1], self.wx[1] \
        , self.nx[0], self.maxnptsx[0]+1 \
        , self.nx[1], self.maxnptsx[1]+1)

    def _save_grid_3D(self):
        self.com.pyfem.pyfem_create_from_tensor_product_3d(self.current_grids \
        , self.id \
        , self.x[0], self.wx[0] \
        , self.x[1], self.wx[1] \
        , self.x[2], self.wx[2] \
        , self.nx[0], self.maxnptsx[0]+1 \
        , self.nx[1], self.maxnptsx[1]+1 \
        , self.nx[2], self.maxnptsx[2]+1)

    def get_real_elts(self):
        # default
        if self.tensorlevel==1:
            return _np.arange(1,self.nel+1)

        if self.tensorlevel==2:
            li_nel = 1

            li_d = 0
            li_n = self.patch.shape[li_d]
            li_p = self.patch.degree[li_d]
            li_un = len(_np.unique(self.patch.knots[li_d]))
            ll_isperiodic = (li_un == li_n + li_p + 1)
            if ll_isperiodic :
                list_u = _np.array([self.patch.knots[li_d][i] for i in range(li_p,li_n+1) ])
                li_nel *= li_n - li_p
            else:
                list_u = _np.unique(self.patch.knots[li_d])
                li_nel *= len(list_u)-1

            for li_d in range (1,self.dim):
                li_n = self.patch.shape[li_d]
                li_p = self.patch.degree[li_d]
                li_un = len(_np.unique(self.patch.knots[li_d]))
                ll_isperiodic = (li_un == li_n + li_p + 1)
                if ll_isperiodic :
                    list_u = _np.array([self.patch.knots[li_d][i] for i in range(li_p,li_n+1) ])
                    li_nel *= li_n - li_p
                else:
                    list_u = _np.unique(self.patch.knots[li_d])
                    li_nel *= len(list_u)-1

            return _np.arange(1,li_nel+1)


    def gen_dbasis(self, dbasis=None, nderiv=1):
        # ...
        # Default dbasis
        # ...
        from scipy.interpolate import splev
        def dxbspline(d,i,x, der):
            """evaluate b-spline first derivative starting at node i at x, d is the direction"""
            lpr_t = np.asarray(self.patch.knots[d])
            li_k = self.patch.degree[d]
            c=np.zeros_like(lpr_t)
#            c=np.zeros(self.patch.shape[d]+self.patch.degree[d]+1)
            if i < 0:
                c[li_n + i]=1.
            else:
                c[i]=1.
            return splev(x,(lpr_t,c,li_k), der=der)
        # ...

        if dbasis is None:
            func = dxbspline
        else:
            func = dbasis

        li_grids_id = self.current_grids
        li_id = self.id

        li_maxnderiv = nderiv
        li_maxp = max(self.patch.degree)
        li_maxni = max(self.nx)
        li_maxdirnpts = self.dirmaxnpts

#        print "====="

#        print li_maxnderiv, li_maxp, li_maxni, li_maxdirnpts
#        print "====="


        # ...
        for li_d in range(0, self.dim):
            li_nderiv = nderiv
            li_p  = self.patch.degree[li_d]
            li_ni = self.nx[li_d]
            li_dirnpts = self.maxnptsx[li_d]

#            print li_nderiv, li_p, li_ni, li_dirnpts

            lpr_dbasisatx = np.zeros((li_nderiv+1, li_p+1, li_dirnpts,li_ni))

            list_x = self.x[li_d]
            ielt = 0
            for x in list_x:
                for i in range(ielt,ielt+li_p+1):
                    for der in range(0,li_nderiv+1):
                        dxbasisatx  = func(li_d,i, x[1:], der)
                        lpr_dbasisatx[der, i-ielt, :, ielt] = dxbasisatx[:]

                ielt += 1

            self.com.pyfem.pyfem_set_tensor_basis ( li_grids_id, li_id, li_d+1 \
            , lpr_dbasisatx \
            , li_nderiv, li_p, li_ni, li_dirnpts \
            , li_maxnderiv, li_maxp, li_maxni, li_maxdirnpts )
        # ...




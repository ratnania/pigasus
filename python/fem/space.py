# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['space']
__date__ ="$Jan 12, 2012 8:56:37 AM$"

# ...
import numpy as _np
from . import common_obj as _com
from . import constants as _cst
from .pigasusObject import *
from .ElementsManager import ElementsManager
from numpy import zeros, asarray
# ...

meshgrid = _np.meshgrid

class space(pigasusObject):
    def __init__(self, as_file=None, geometry=None \
    , mapping = None, tensor=True, tensorlevel=1, al_storeddata=False):
        pigasusObject.__init__(self)

        self.type="space"
        self.tensor=tensor
        self.tensorlevel=tensorlevel
        self.grids = None
        self.mapping = None
        self.metric = None

#        # needed for assembly, when evaluating fields
#        self._currentPatch = None
#        self._currentPatchID = None

        self.list_faces_dir = []
        self.list_faces_duplicated = []
        self.list_faces_duplicata  = []

        # physical points
        self._points = None
        # parametric (sites) points
        self._sites  = None

        if as_file is not None:
            import caid.cad_geometry as cg
            self.geometry = cg.cad_geometry(file=as_file)
        else :
            if geometry is None:
                print("SERIOUS ERROR while constructing the space")
                import sys
                sys.exit(2)
            self.geometry = geometry

        self.EM = ElementsManager(self.geometry)

        self.dim = self.geometry.dim
        self.Rd = self.geometry.Rd
        self.ndof = 1
        self.storeddata = al_storeddata
        self.nderiv_pts = 1 # used to compute the position, the 1st derivative

        self.maxnen = 0

        self.connectivity = None
        self.boundary_conditions = None

        self.exterior_map = False
        if mapping is not None:
            self.exterior_map = True

        self.usemetric = False
        self.metric_id = -1
        self.metric = None

#        self.updatePoints = True

        from caid.numbering.connectivity import connectivity
        from caid.numbering.boundary_conditions import boundary_conditions
        self.connectivity = connectivity(self.geometry, ai_ndof=self.ndof)
        self.boundary_conditions = boundary_conditions(self.geometry)

        if self.dim not in [1,2,3] :
            print("Error space: dimension not implemented yet")

        # setting the maxnen for the space
        self.maxnen = max (self.connectivity.list_nen)

        # this must be the last thing to do
        self.id = self.com.nspaces
        self.com.nspaces += 1
        self.com.spaces.append(self)

    def set_currentPatch(self, patch):
        self.grids.set_currentPatch(patch)

    @property
    def currentPatch(self):
        return self.grids.currentPatch

    @property
    def currentPatchID(self):
        return self.grids.currentPatchID

    def setInfoData(self):
        """
        prints informations about the current matrix
        """
        self.infoData['id'] = str(self.id)
        self.infoData['size'] = str(self.size)
        self.infoData['tensor'] = str(self.tensor)
        self.infoData['tensorlevel'] = str(self.tensorlevel)
        self.infoData['storeddata'] = str(self.storeddata)
        self.infoData['exterior_map'] = str(self.exterior_map)
        self.infoData['dim'] = str(self.dim)
        self.infoData['ndof'] = str(self.ndof)
        self.infoData['maxnen'] = str(self.maxnen)

    def set_points(self, values):
        self._points = values

    def get_points(self):
        return self._points

    def set_sites(self, values):
        self._sites = values

    def get_sites(self):
        return self._sites

    def add_patchs(self):
        # ...
        from caid.io import formatter
        fmt = formatter()
        def _add_patch(ai_id, nrb):
            """
            this routine updates the current patch to
            the geometries structure defined in Fortran
            """
            li_dim = nrb.dim
            lpi_N = nrb.shape
            lpi_P = nrb.degree
            li_maxnk = max(lpi_N)+max(lpi_P)+1

            lpr_u = _np.zeros((li_maxnk, li_dim), dtype=_np.double)

            for d in range(0,li_dim):
                li_N = lpi_N[d]
                li_P = lpi_P[d]
                lpr_u [0:li_N+li_P+1, d] = nrb.knots[d]

            lpr_P = fmt.to_list(nrb.points , nrb.dim, nrb.shape, Rd=self.Rd)
            lpr_W = fmt.to_list(nrb.weights, nrb.dim, nrb.shape, Rd=1)

            if nrb.rational:
                li_rational = 1
            else:
                li_rational = 0

            self.com = _com.common_obj()
            self.com.pyfem.set_patch_space(self.id, ai_id, li_rational, lpi_N, lpi_P, lpr_u, lpr_P, lpr_W)
        # ...
        # ...
        def _add_patchs():
            li_nmp = self.geometry.npatchs

            for li_id in range(0,li_nmp):
                lo_patch = self.geometry[li_id]
    #            lo_patch.print_patch()
                _add_patch(li_id, lo_patch)
        # ...

        li_nmp = self.geometry.npatchs
        for li_id in range(0,li_nmp):
            lo_patch = self.geometry[li_id]
#            lo_patch.print_patch()
            _add_patch(li_id, lo_patch)

    def add_connectivity(self):
        npatchs = self.grids.npatchs
        maxnen  = self.maxnen
        maxnelt = max ([self.grids.list_grid[i].nel for i in range(0, npatchs)])
        nrealelt = max([self.connectivity.LM[i].shape[1] for i in range(0, npatchs)])

        RealElts = _np.zeros((npatchs, maxnelt) ,dtype=_np.int)
#        print "======="
#        print self.id
#        print "maxnelt ", maxnelt
        LM  = _np.zeros((npatchs, maxnen, nrealelt) ,dtype=_np.int)
        IEN = _np.zeros((npatchs, maxnen, nrealelt) ,dtype=_np.int)

        # ...
        _RealElts = self.grids.get_real_elts()
#        for i in range(0, len(_RealElts)):
#            print _RealElts[i].shape
        for i in range(0, npatchs):
            _nelt = len(_RealElts[i])
            RealElts[i,0:_nelt] = _RealElts[i][0:_nelt]
#            for j in range(0, maxnelt):
#                RealElts[i,j] = _RealElts[i][j]
        # ...

        # ...
        for i in range(0, npatchs):
            _maxnen, _nrealelt = self.connectivity.LM[i].shape
            LM[i,0:_maxnen,0:_nrealelt]  = self.connectivity.LM[i][0:_maxnen,0:_nrealelt]
            IEN[i,0:_maxnen,0:_nrealelt] = self.connectivity.IEN[i][0:_maxnen,0:_nrealelt]
#            for j in range(0, _maxnen):
#                for k in range(0, _nrealelt):
#                    LM[i,j,k]  = self.connectivity.LM[i][j,k]
#                    IEN[i,j,k] = self.connectivity.IEN[i][j,k]
        # ...

        ID = self.connectivity.ID

        self.com.pyfem.pyfem_create_connectivity(self.id, npatchs, maxnen, maxnelt)
        self.com.pyfem.pyfem_set_connectivity_id(self.id, ID)
        self.com.pyfem.pyfem_set_connectivity_real_elts(self.id, RealElts)
        self.com.pyfem.pyfem_set_connectivity_lm(self.id, LM)
        self.com.pyfem.pyfem_set_connectivity_ien(self.id, IEN)

#        print ID
#        print IEN
#        print LM
#        print RealElts

        for li_patch in range(0, npatchs):
            ID_loc = self.connectivity.ID_loc[li_patch]
#            print "ID_loc : ", ID_loc
            self.com.pyfem.pyfem_set_connectivity_loc(self.id, li_patch, ID_loc)

#        , self.connectivity.IEN, self.connectivity.LM)
#        print "self.connectivity.ID=", self.connectivity.ID
#        print "self.connectivity.IEN=", self.connectivity.IEN
#        print "self.connectivity.LM=", self.connectivity.LM
#        print "self.grids.get_real_elts()=",self.grids.get_real_elts()

    def create_grids(self, k=None, type=None, space=None \
    , profile="volume", metric = None \
    , list_nodes=None, list_w=None, list_faces=None ):
        if space is None:
            from . import grids as gr
            if list_nodes is None:
                if metric is not None:
                    self.metric_id = metric.id
                    self.usemetric = True
                    self.metric = metric
                    if metric.type == "cad_mapping":
                        self.exterior_map = True

                self.grids = gr.grids(profile, api_k=k, ao_geometry=self.geometry \
                , tensorlevel=self.tensorlevel, as_type=type, ai_ndof=self.ndof,
                                      metric = metric, list_faces=list_faces)
            else:
                self.grids = gr.grids(profile, api_k=k, ao_geometry=self.geometry \
                , tensorlevel=self.tensorlevel, as_type=type, ai_ndof=self.ndof \
                , list_nodes = list_nodes, list_w = list_w, list_faces=list_faces)
        else :
            if self.ndof > space.grids.ndof:
                space.grids.ndof = self.ndof
            self.grids = space.grids

    def get_points_advanced(self, ai_patch_id):
        V       = self
        geo     = V.geometry
        G       = V.grids
        nderiv  = V.nderiv_pts
#        print "Current space ", self.id
#        print "nderiv = ", nderiv

        patch_id = ai_patch_id
        zeros = _np.zeros

        _sites = G.get_sites()[patch_id]

        npts = G.k

        srf = geo[patch_id]
        nx = G.list_grid[patch_id].nx

        # ... defining the evaluate function
#        EVALUATE = srf.evaluate_deriv
        def EVALUATE(u=None, v=None, w=None):
            return srf.evaluate_deriv(u=u, v=v, w=w, nderiv=nderiv)
#        APPEND = True
        APPEND = False
        # ...

        from .utils import evaluator
        X = evaluator(srf.dim, _sites, EVALUATE, APPEND, nx=nx, npts=npts,
                      fmt=True)
        return X

    def computePoints(self, patch_id):
#        print "Enter computePoints"
        li_patch = patch_id
        lpr_parampts = self.get_parametricPoints(ai_patch_id=li_patch)

        _G = self.grids.list_grid[li_patch]
        shape_pts = [_G.nel, _G.npts]

        if self.metric is not None:
            if self.metric.type == "analytic":
                lpr_pts = self.metric.update_analytic(li_patch, lpr_parampts, _G.nel, _G.npts)
            # ------------
            if self.metric.type == "cad_geometry":
                V       = self
                G       = V.grids

                nderiv  = V.nderiv_pts
                sites = G.get_sites()[li_patch]
                k = G.k
                nx = G.list_grid[li_patch].nx
                nel = _G.nel
                npts = _G.npts

                lpr_pts = self.metric.update_cad_geometry(li_patch, sites, nderiv, k, nx, nel, npts)

        else:
            lpr_pts  = self.get_points_advanced(ai_patch_id=li_patch)
#        print "******"
#        print lpr_pts.shape

        self.set_sites(lpr_parampts)
        self.set_points(lpr_pts)

    def get_parametricPoints(self, ai_patch_id):
        G = self.grids.list_grid[ai_patch_id]
        sts = self.com.pyfem.get_space_parametricpoints(self.id, ai_patch_id    \
    , G.nel, G.maxnpts, G.dim )

        n = sts.shape
        list_sts = []
        for k in range(0,n[2]):
            new_sts = _np.zeros(n[0]*n[1])
            for i in range(0,n[0]):
                for j in range(0,n[1]):
                    new_sts[n[1]*i+j] = sts[i,j,k]
            list_sts.append(new_sts)

        return list_sts

    def get_weightsGrids(self, ai_patch_id):
        li_nel = self.grids.list_grid[ai_patch_id].nel
        li_npts = self.grids.list_grid[ai_patch_id].maxnpts
        return self.com.pyfem.pyfem_get_weights_grid(self.grids.id, ai_patch_id, li_nel, li_npts)

    def dirichlet(self, faces):
        if len(faces) > 0 :
            self.boundary_conditions.dirichlet(self.geometry, faces)
            self.list_faces_dir = faces

    def duplicate(self, \
                  faces_base=None, faces=None, \
                  isPeriodic=None):
        if isPeriodic is None:
            pass # TODO
        _faces_base = []
        _faces = []
        for (pf_base, pf, periodic) in zip(faces_base, faces, isPeriodic):
#            print "pf_base =", pf_base, "       pf =", pf
            pbase_id = pf_base[0] ; fbase_id = pf_base[1]
            p_id = pf[0] ; f_id = pf[1]
            if periodic:
#                print ">>> is Periodic"
                nrb_base = self.geometry[pbase_id]
                nrb      = self.geometry[p_id]
                assert(nrb_base.dim == nrb.dim)
                axis = None
                if nrb_base.dim == 1:
                    if f_id in [0,1]:
                        axis = 0
                if nrb_base.dim == 2:
                    if f_id in [0,2]:
                        axis = 0
                    if f_id in [1,3]:
                        axis = 1
                if nrb_base.dim == 3:
                    # TODO
                    print ("Not yet implemented")

                degree = nrb_base.degree[axis]
                for i_base in range(0, degree):
                    i = i_base - degree + 1
                    new_pf_base = [pbase_id, fbase_id, i_base]
                    new_pf      = [p_id, f_id, i]
                    _faces_base.append(new_pf_base)
                    _faces.append(new_pf)
            else:
                new_pf_base = [pbase_id, fbase_id, 0]
                new_pf      = [p_id, f_id, 0]
                _faces_base.append(new_pf_base)
                _faces.append(new_pf)

        self.boundary_conditions.duplicate(self.geometry)
        self.list_faces_duplicated = _faces_base
        self.list_faces_duplicata  = _faces

    def set_boundary_conditions(self):
        self.connectivity.init_data_structure(self.boundary_conditions)
        self.size = self.connectivity.size

    def get_dir_nel(self):
        """
        returns the number of elements for each direction
        """
        import numpy as np
        list_nel = []
        li_npatchs = self.grids.npatchs
        for li_id in range(0,li_npatchs):
            list_localnel = []
            lo_grid = self.grids.list_grid[li_id]
            if lo_grid.dim == 2:
                list_localnel = [lo_grid.nx, lo_grid.ny]
            else:
                print("get_dir_nel : Not yet implemented")
            list_nel.append(_np.array(list_localnel))
        return _np.array(list_nel)

    def to_meshGridFormat(self, ai_patch, apr_values):
        nel = self.grids.list_grid[ai_patch].nel
        li_dim = self.dim
        if li_dim not in [1,2] :
            print("to_meshGridFormat not yet implemented for dimension = ", li_dim)
        if li_dim == 1 :
            return self._to_meshGridFormat_1D(ai_patch, apr_values)
        if li_dim == 2 :
            return self._to_meshGridFormat_2D(ai_patch, apr_values)

    def _to_meshGridFormat_1D(self, ai_patch, apr_values):

        li_ntotal = 0
        list_npts = []
        for x in self.grids.list_grid[ai_patch].x[0]:
            list_npts.append(int(x[0]))
            li_ntotal += _np.asarray(x[0])

        li_nf = len(apr_values)
#        li_nf = apr_values.shape[0]
        F = _np.zeros((li_nf, li_ntotal))

        li_nel = len(_np.unique(_np.asarray(self.geometry[ai_patch].knots[0]))) - 1

        nel = li_nel
        npts = self.grids.list_grid[ai_patch].maxnpts
        list_val = []
        for lv in apr_values:
            list_val.append(lv.reshape([nel,npts]))
        elt = 0
        istart = 0
        for I in range(0, li_nel):
            nx = list_npts[I]
            for d in range(0,li_nf):
                F[d,istart:istart+nx] = list_val[d][elt,:]
            istart += nx
            elt += 1

        return F

    def _to_meshGridFormat_2D(self, ai_patch, apr_values):
        # first we compute the N1,N2 total points in each direction
        lpi_ntotal = _np.asarray([0,0])
        list_npts = []
        i = 0
        for Lx in self.grids.list_grid[ai_patch].x:
            list_nx = []
            for x in Lx:
                list_nx.append(int(x[0]))
                lpi_ntotal[i] += _np.asarray(x[0])
            list_npts.append(list_nx)
            i += 1

        li_nf = len(apr_values)
#        li_nf = apr_values.shape[0]
        F = _np.zeros((li_nf, lpi_ntotal[0], lpi_ntotal[1]))

        lpi_nel = [0,0]
        lpi_nel[0] = len(_np.unique(_np.asarray(self.geometry[ai_patch].knots[0]))) - 1
        lpi_nel[1] = len(_np.unique(_np.asarray(self.geometry[ai_patch].knots[1]))) - 1

        nel  = lpi_nel[0] * lpi_nel[1]
        npts = self.grids.list_grid[ai_patch].maxnpts
        list_val = []
        for lv in apr_values:
            list_val.append(lv.reshape([nel,npts]))

        for J in range(0, lpi_nel[1]):
            ny = list_npts[1][J]
            jstart = sum(_np.asarray(list_npts[1][0:J]))
            for I in range(0, lpi_nel[0]):
                istart = sum(_np.asarray(list_npts[0][0:I]))
                elt = J * lpi_nel[0] + I
                nx = list_npts[0][I]
                lpi_npts = [nx,ny]
                for d in range(0,li_nf):
                    lpr_val = list_val[d][elt,:].reshape(lpi_npts[::-1]).transpose()
                    F[d,istart:istart+nx,jstart:jstart+ny] = lpr_val[0:nx,0:ny]
        list_F = []
        for d in range(0, li_nf):
            list_F.append(F[d,:,:].transpose())
        return list_F


    def set_weights(self, type=None):

        li_patch_id = 0
        for geo in self.geometry:
            list_t = []
            for li_d in range(0,self.dim):
                if type[li_d] is None :
                    ls_curtype = "N"
                else :
                    ls_curtype = type[li_d]
                if ls_curtype=="N":
                    li_n = geo.shape[li_d]
                    lpr_t = _np.ones(li_n,dtype=_np.double)
                if ls_curtype=="T":
                    li_n = geo.shape[li_d]
                    lpr_t = _np.ones(li_n,dtype=_np.double)
                    lpr_t *= 2.0
                if ls_curtype=="D":
                    li_n = geo.shape[li_d]
                    li_k = geo.degree[li_d] + 1
                    lpr_knots = _np.array(geo.knots[li_d])
                    lpr_t = li_k / (lpr_knots[0+li_k:li_n+li_k] - lpr_knots[0:li_n])

                self.com.pyfem.set_space_weights(self.id,li_patch_id, li_d+1 \
                , lpr_t, li_n)
                print("lpr_t=", lpr_t)
#                list_t.append(lpr_t)
#            if self.dim == 1 :
#                lpr_weights = _np.array(list_t[0])
#            if self.dim == 2 :
#                lpr_a = _np.array(list_t[0])
#                lpr_b = _np.array(list_t[1])
#                lpr_weights = _np.zeros(len(lpr_a)*len(lpr_b), dtype=_np.double)
#                li_index = 0
#                for b in lpr_b:
#                    for a in lpr_a:
#                        lpr_weights[li_index] = a * b
#                        li_index += 1
#            if self.dim == 3 :
#                lpr_a = _np.array(list_t[0])
#                lpr_b = _np.array(list_t[1])
#                lpr_c = _np.array(list_t[2])
#                lpr_weights = _np.zeros(len(lpr_a)*len(lpr_b)*len(lpr_c), dtype=_np.double)
#                li_index = 0
#                for c in lpr_c:
#                    for b in lpr_b:
#                        for a in lpr_a:
#                            lpr_weights[li_index] = a * b * c
#                            li_index += 1
#            print "lpr_weights=", lpr_weights
#            self.com.pyfem.set_space_weights(self.id,li_patch_id \
#            , lpr_weights, li_index)
            li_patch_id += 1

    def update_geometry(self, geometry):
        # ...
        from caid.io import formatter
        fmt = formatter()
        def update_points_patch(ai_id, nrb):
            """
            this routine updates the current patch to
            the geometries structure defined in Fortran
            """
            li_dim = nrb.dim
            lpi_N = nrb.shape
            lpi_P = nrb.degree
            li_maxnk = max(lpi_N)+max(lpi_P)+1

            lpr_u = _np.zeros((li_maxnk, li_dim), dtype=_np.double)

            for d in range(0,li_dim):
                li_N = lpi_N[d]
                li_P = lpi_P[d]
                lpr_u [0:li_N+li_P+1, d] = nrb.knots[d]

            lpr_P = fmt.to_list(nrb.points , nrb.dim, nrb.shape, Rd=self.Rd)
            lpr_W = fmt.to_list(nrb.weights, nrb.dim, nrb.shape, Rd=1)

            self.com = _com.common_obj()
            self.com.pyfem.update_points_patch_space(self.id ,ai_id, lpr_P, lpr_W)
        # ...
        self.geometry = geometry
        li_nmp = self.geometry.npatchs
        for li_id in range(0,li_nmp):
            lo_patch = geometry[li_id]
            update_points_patch(li_id, lo_patch)

    def inv_vect(self, ai_patch, apr_values):
        """
        this routine returns the input vector apr_values into a matrix
        apr_values is a vector on the space V
        the output is a  2D or 3D matrix
        """

        lpi_nV = self.geometry[ai_patch].shape
        li_dim = self.dim
        lpi_istart  = [0]*li_dim
        li_Astart   = self.connectivity.get_A_ind(ai_patch, lpi_istart)
        lpi_iend    = _np.asarray(lpi_nV)-1
        li_Aend     = self.connectivity.get_A_ind(ai_patch, lpi_iend)

        list_val = []
        for li_A in range(li_Astart, li_Aend + 1):
            li_P = self.connectivity.ID[li_A]
            if li_P == 0:
                list_val.append(0.0)
            else:
                list_val.append(apr_values[li_P-1])

        lpr_val = _np.array(list_val).reshape(lpi_nV[::-1]).transpose()

        return lpr_val

    def patch_coef(self, ai_patch, apr_values):
        """
      this routine gives the coefs for the given patch, from the global coef vector
        """

        lpi_nV = self.geometry[ai_patch].shape
        li_dim = self.dim
        lpi_istart  = [0]*li_dim
        li_Astart   = self.connectivity.get_A_ind(ai_patch, lpi_istart)
        lpi_iend    = _np.asarray(lpi_nV)-1
        li_Aend     = self.connectivity.get_A_ind(ai_patch, lpi_iend)

        list_val = []
        for li_A in range(li_Astart, li_Aend + 1):
            li_P = self.connectivity.ID[li_A]
            if li_P == 0:
                list_val.append(0.0)
            else:
                list_val.append(apr_values[li_P-1])

        lpr_val = _np.array(list_val)

        return lpr_val




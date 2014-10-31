# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['grids']
__date__ ="$Jan 13, 2012 4:22:57 PM$"

from . import common_obj as _com
from .pigasusObject import *

class grids(pigasusObject):
    def __init__(self, profile, api_k=None, ao_geometry=None    \
    , tensor=True, tensorlevel=1, as_type="lobatto", ai_ndof=1 \
    , metric = None, list_nodes=None, list_w=None, list_faces=None ):
        pigasusObject.__init__(self)

        if api_k is None:
            print("grids : api_k must be given")
            import sys
            sys.exit(2)

        if ao_geometry is None:
            print("grids : SERIOUS ERROR while constructing the grid")
            import sys
            sys.exit(2)

        self.geometry = ao_geometry
        self.k = [k+1 for k in api_k]

        self.list_grid = []
        self.tensorlevel = tensorlevel

        self.usemetric = False
        self.metric_id = 0
        self.metric = None
        if metric is not None:
            self.metric_id = metric.id
            self.usemetric = True
            self.metric = metric

        self.npatchs = self.geometry.npatchs

        self.id = self.com.ngrids
        for li_id in range(0,self.geometry.npatchs):
            lo_patch = self.geometry[li_id]
            if profile=="volume":
                from . import grid as gr
                if list_nodes is None:
                    lo_grid = gr.grid(self.id, li_id, lo_patch, api_k, tensorlevel, as_type)
                else:
                    lo_grid = gr.grid(self.id, li_id, lo_patch, api_k, tensorlevel, as_type \
                    , list_nodes = list_nodes[li_id], list_w = list_w[li_id])
            if profile=="boundary":
                from . import boundary_grid as gr
                lo_grid = gr.boundary_grid(self.id, li_id, lo_patch, api_k,
                                           tensor, as_type, faces=list_faces[li_id])
            if profile=="surface":
                print("surface grid not yet validated")
                import sys; sys.exit(0)
            if profile=="vertex":
                print("vertex grid not yet validated")
                import sys; sys.exit(0)
            self.list_grid.append(lo_grid)

        # DOMAIN'S DIMENSION
        self.dim = self.list_grid[0].dim
        self.Rd = self.list_grid[0].Rd
        self.ndof = ai_ndof
        # for each patch we store the max-number of points in each element
        self.maxnpts = max( [G.maxnpts for G in self.list_grid])
        # for each patch we store the max-number of points in each element among
        # all directions
        self.dirmaxnpts = max( [G.dirmaxnpts for G in self.list_grid])
        # for each patch we store if we can use a tensor product for assembling basis
        self.tensor = tensor

        self._patchs_toassembly = None
        self._elts_toassembly = None

        # needed for assembly, when evaluating fields
        self._currentPatch = None
        self._currentPatchID = None

        # this must be the last thing to do
        self.com.ngrids += 1
        self.com.grids.append(self)
        self.com.list_Grfields_id.append([])
        self.com.list_Groperators_id.append([])
        self.com.list_Grnorms_id.append([])

    def set_currentPatch(self, patch):
        self._currentPatch = patch
        self._currentPatchID = self.geometry.index(patch)

    @property
    def currentPatch(self):
        return self._currentPatch

    @property
    def currentPatchID(self):
        return self._currentPatchID

    def get_sites(self, patch_id=None):
        if patch_id is None:
            list_patchs = list(range(0, self.npatchs))
        else:
            list_patchs = [patch_id]
        list_x = []
        for id in list_patchs:
            list_x.append(self.list_grid[id].grid1D)
        return list_x

    def get_fields_id(self):
        return self.com.list_Grfields_id[self.id]

    def add_field_id(self, F):
        self.com.list_Grfields_id[self.id].append(F.id)
        # return the local id of the field
        return len(self.com.list_Grfields_id[self.id]) - 1

    def get_operators_id(self):
        return self.com.list_Groperators_id[self.id]

    def add_operator_id(self, M):
        self.com.list_Groperators_id[self.id].append(M.id)
        # return the local id of the matrix
        return len(self.com.list_Groperators_id[self.id]) - 1

    def get_norms_id(self):
        return self.com.list_Grnorms_id[self.id]

    def add_norm_id(self, N):
        self.com.list_Grnorms_id[self.id].append(N.id)
        # return the local id of the norm
        return len(self.com.list_Grnorms_id[self.id]) - 1

    def print_grids(self, ai_id=None):
        if ai_id is not None:
            print("******")
            print("Grid id =", ai_id)
            print("******")
            self.com.pyfem.pyfem_print_grid(self.id, ai_id)
        else:
            for li_id in range(0,self.npatchs):
                print("******")
                print("Grid id =", li_id)
                print("******")
                self.com.pyfem.pyfem_print_grid(self.id, li_id)

    def save_grids(self):
         for li_id in range(0,self.npatchs):
            lo_grid = self.list_grid[li_id]
            lo_grid.save_grid()

    def get_real_elts(self):
        list_elts = []
        for li_id in range(0,self.npatchs):
            lo_grid = self.list_grid[li_id]
            list_elts.append(lo_grid.get_real_elts())
        return list_elts

    def gen_dbasis(self, dbasis=None, nderiv=1):
        for G in self.list_grid:
            G.gen_dbasis(dbasis=dbasis, nderiv=nderiv)

    def set_patchs_toassembly(self, list_patchs):
        if list_patchs is None:
            self._patchs_toassembly = list(range(0, self.npatchs))
        else:
            self._patchs_toassembly = list_patchs

    def get_patchs_toassembly(self):
        return self._patchs_toassembly

    def set_elts_toassembly(self, list_elts):
        if list_elts is None:
            self._elts_toassembly = [] #this means that we will assemble all elements of the current patch
        else:
            self._elts_toassembly = list_elts

    def get_elts_toassembly(self):
        return self._elts_toassembly

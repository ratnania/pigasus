# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__date__ ="$Jan 12, 2012 9:32:06 AM$"
__all__ = "mapping"

from . import common_obj as _com
import numpy as _np
from .common_obj import *
from .pigasusObject import *

class mapping(pigasusObject):
    def __init__(self, geometry, tensor=True, al_storeddata=False):
        lo_com = common_obj()

        self.type="mapping"
        self.tensor=tensor
        self.storeddata = al_storeddata

        self.geometry = geometry

        # this must be the last thing to do
        self.id = lo_com.nmappings
        lo_com.nmappings += 1
        lo_com.mappings.append(self)

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

            lpr_P = fmt.to_list(nrb.points , nrb.dim, nrb.shape, Rd=3)
            lpr_W = fmt.to_list(nrb.weights, nrb.dim, nrb.shape, Rd=1)

            if nrb.rational:
                li_rational = 1
            else:
                li_rational = 0

            self.com = _com.common_obj()
            self.com.pyfem.set_patch_mapping(self.id ,ai_id, li_rational, lpi_N, lpi_P, lpr_u, lpr_P, lpr_W)
        # ...
        # ...
        def _add_patchs():
            li_nmp = self.geometry.npatchs

            for li_id in range(0,li_nmp):
                lo_patch = self.geometry[li_id]
                _add_patch(li_id, lo_patch)
        # ...

        li_nmp = self.geometry.npatchs
        for li_id in range(0,li_nmp):
            lo_patch = self.geometry[li_id]
            _add_patch(li_id, lo_patch)

# -*- coding: UTF-8 -*-
#! /usr/bin/python

__author__="ARA"
__all__ = ['pastix']
__date__ ="$Jan 4, 2013 3:12:55 PM$"

from . import common_obj as _com
from . import constants as _cst
import numpy as _np
from .pigasusObject import *

class pastix(pigasusObject):
    def __init__(self):
        pigasusObject.__init__(self)
        print("__init__ : TODO")

    def set_rhs(self, matrix, field=None, apr_val=None):
        if field is not None:
            self._set_rhs_field(matrix, field)
        if apr_val is not None:
            self._set_rhs_array(matrix, apr_val)

    def get_solution(self, matrix, field=None):
        if field is not None:
            self._get_solution_field(matrix, field)
        else:
            return self._get_solution_array(matrix)

    def _set_rhs_field(self, matrix, field):
        print("_set_rhs_id : TODO")
        self.com.pyfem.solver_setfieldrhs(matrix.id, field.id )

    def _set_rhs_array(self, matrix, apr_val):
        print("_set_rhs_array : TODO")
        self.com.pyfem.solver_setvaluesrhs(matrix.id, apr_val )

    def _get_solution_field(self, matrix, field):
        print("_get_solution_id : TODO")
        li_size = matrix.spaces[0].size
        self.com.pyfem.solver_getfieldsolution ( matrix.id, field.id )

    def _get_solution_array(self, matrix):
        print("_get_solution_array : TODO")
        li_size = matrix.spaces[0].size
        return self.com.pyfem.solver_getsolution ( matrix.id , li_size )

# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="root"
__date__ ="$Jan 12, 2012 2:26:57 PM$"

# ***************************************************
"""
SOLVERS
this is a disctionary where for each Library we give the corresponding solvers
"""
# ***************************************************
# ... solvers
SPM_SOLVER_CUSP_CG                    = 0
SPM_SOLVER_CUSP_GMRES                 = 1
SPM_SOLVER_CUSP_BICGSTAB              = 2
SPM_SOLVER_CUSP_JACOBI                = 3
SPM_SOLVER_BASIC_CG                   = 4
SPM_SOLVER_BASIC_GS                   = 5
# ***************************************************

# ***************************************************
"""
TYPE MATRICES
"""
# ***************************************************
GENERIC_MATRIX      = 0
BLOCK_MATRIX        = 1
# ***************************************************

# ***************************************************
"""
TYPE OPERATORS
"""
# ***************************************************
MASS                = 0
STIFFNESS           = 1
ADVECTION           = 2
SECOND_DERIV        = 3
# ***************************************************

# ***************************************************
"""
TYPE PROJECTORS
"""
# ***************************************************
PROJECTION_L2   = 1
FIELD_OPERATOR  = 2
# ***************************************************

# ***************************************************
"""
TYPE DIAGNOSTICS
"""
# ***************************************************
IDENTITY         = 1
DIAGS_ADVECTION  = 2
# ***************************************************

# ***************************************************
"""
TYPE NORMS
"""
# ***************************************************
NORM_L2 = 1
NORM_H1 = 2
# ***************************************************

# ***************************************************
"""
FIELDS OPERATORS TYPES
"""
# ***************************************************
IDENTITY               = 1
GRAD                   = 2
CURL                   = 3
IDENTITY_VECT          = 4
SECOND_DERIV_FIELD     = 5
GRAD_S                 = 6
SECOND_DERIV_S_FIELD   = 8
HESSIAN_FIELD          = 9

list_FIELD_OPERATORS = [GRAD, CURL, SECOND_DERIV_FIELD \
, GRAD_S, SECOND_DERIV_S_FIELD, HESSIAN_FIELD]

def MAP_FIELD_DIMOUT(dim, ID):
    if dim ==1:
        if ID==IDENTITY:
            return 1
        if ID==IDENTITY_VECT:
            return 1
        if ID==GRAD_S:
            return 1
        if ID==GRAD:
            return 1
        if ID==CURL:
            return 1
        if ID==SECOND_DERIV_FIELD:
            return 1
        if ID==SECOND_DERIV_S_FIELD:
            return 1
        if ID==HESSIAN_FIELD:
            return 1
    if dim ==2:
        if ID==IDENTITY:
            return 1
        if ID==IDENTITY_VECT:
            return 1
        if ID==GRAD_S:
            return 1
        if ID==GRAD:
            return 2
        if ID==CURL:
            return 2
        if ID==SECOND_DERIV_FIELD:
            return 3
        if ID==SECOND_DERIV_S_FIELD:
            return 1
        if ID==HESSIAN_FIELD:
            return 1
# ***************************************************

# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__date__ ="$Jul 27, 2012 11:35:03 PM$"

import numpy as _np
def internal_func (func, ai_dim):
    """
    this routine return a function-type that is used in field/matrix/...
    so that the user can call it with only a list of points (given by get_points)
    """
    def lo_func_1D (P) :
        x = _np.asarray(P[0][:])
        return _np.asarray([x[:] - x[:] + z for z in _np.asarray(func(x))])

    def lo_func_2D (P) :
        x = _np.asarray(P[0][:])
        y = _np.asarray(P[1][:])
        return _np.asarray([x[:] - x[:] + z for z in _np.asarray(func(x,y))])

    def lo_func_3D (P) :
        x = _np.asarray(P[0][:])
        y = _np.asarray(P[1][:])
        z = _np.asarray(P[2][:])
        return _np.asarray([x[:] - x[:] + z for z in _np.asarray(func(x,y,z))])

    if ai_dim == 1:
        return lo_func_1D
    elif ai_dim == 2:
        return lo_func_2D
    elif ai_dim == 3:
        return lo_func_3D
    else :
        print("func_tools : Not yet implemented")
        import sys; sys.exit(0)

#def internal_arg_func (func, ai_dim):
#    """
#    same as internal_func, but allows an aditional argument, which will be a list of fields
#    todo : pas teste sur matrice stiffness (ou advection, avec dim param >=2)
#    """
#    def lo_func_1D (F,P) :
#        x = _np.asarray(P[0][:])
#        return _np.asarray([x[:] - x[:] + z for z in _np.asarray(func(F,x))])
#
#    def lo_func_2D (F,P) :
#        x = _np.asarray(P[0][:])
#        y = _np.asarray(P[1][:])
#        return _np.asarray([x[:] - x[:] + z for z in _np.asarray(func(F,x,y))])
#
#    def lo_func_3D (F,P) :
#        x = _np.asarray(P[0][:])
#        y = _np.asarray(P[1][:])
#        z = _np.asarray(P[2][:])
#        return _np.asarray([x[:] - x[:] + z for z in _np.asarray(func(F,x,y,z))])
#
#    if ai_dim == 1:
#        return lo_func_1D
#    elif ai_dim == 2:
#        return lo_func_2D
#    elif ai_dim == 3:
#        return lo_func_3D
#    else :
#        print "func_tools : Not yet implemented"
#        import sys; sys.exit(0)

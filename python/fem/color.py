# -*- coding: UTF-8 -*-
#! /usr/bin/python

__author__="ARA"
__all__ = ['color_operator', 'color_field', 'color', 'manager',
           'manager_operators', 'manager_fields']
__date__ ="$Mai 08, 2014 10:50:00 PM$"

from numpy import asarray


class myList:
    def __init__(self):
        self._list = []
        self._currentElt     = -1

    @property
    def list(self):
        return self._list

    def index(self, obj):
        return self._list.index(obj)

    @property
    def n(self):
        return len(self._list)

    def __next__(self):
        if self.n == 0:
            raise StopIteration
        self._currentElt += 1
        if self._currentElt >= self.n:
            self._currentElt = -1
            raise StopIteration
        return self._list[self._currentElt]

    def __iter__(self):
        return self

    def __getitem__(self, key):
        return self._list[key]

    def append(self, obj):
        self._list.append(obj)

    def reset(self):
        self._list = []

from .pigasusObject import *
class color(myList, pigasusObject):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        obj.id = None
        return obj

    def __init__(self, objects=[]):
        pigasusObject.__init__(self)
        myList.__init__(self)

        for obj in objects:
            self.append(obj)
#        print "a new color has been created"

    def append(self, obj):
        if self.is_compatible(obj):
            myList.append(self,obj)
        else:
            print("Warning: Incompatible object with the current color: will not append the current object.")

    def is_compatible(self, obj):
        """
        must be redefined for each color-object
        """
        return True

    def update(self, my_type, my_subtype=-1):
        self.com.pyfem.set_color(self.id, my_type, self.n, my_subtype)
        list_id = asarray([Obj.id for Obj in self._list])
        self.com.pyfem.set_color_objects (self.id, list_id, len(list_id))

class color_operator(color):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        obj.id = None
        return obj

    def __init__(self, operators=[]):
        self._type = None
        self._spaces = None
        color.__init__(self, objects=operators)

        self.id = self.com.ncolors
        self.com.ncolors += 1
        self.com.colors.append(self)
        print("### new color_operator \n")

    def is_compatible(self, oper):
        if oper is None:
            return False
        if self._type is None:
            self._type = oper.type
            self._spaces = oper.spaces
        if self._type == oper.type \
           and self._spaces[0] == oper.spaces[0] \
           and self._spaces[1] == oper.spaces[1] :
            return True
        else:
            return False

    def update(self):
        from .constants import COLOR_OPERATOR
        color.update(self, COLOR_OPERATOR)

class color_field(color):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        obj.id = None
        return obj

    def __init__(self, fields=[]):
        self._type  = None
        self._space = None
        color.__init__(self, objects=fields)

        self.id = self.com.ncolors
        self.com.ncolors += 1
        self.com.colors.append(self)
        print("### new color_field \n")

    def is_compatible(self, F):
        if F is None:
            return False
        if self._type is None:
            self._type  = F.type
            self._space = F.space
        if self._type == F.type \
           and self._space == F.space:
            return True
        else:
            return False

    def update(self):
        from .constants import COLOR_FIELD
        if self._type is not None:
            color.update(self, COLOR_FIELD, self._type)
        else:
            raise("Found a color_field with subtype equal to None.")

from . import common_obj as com
class manager(color):
    def __init__(self, colors=[]):
        self.com = com.common_obj()
        self._type = None
        color.__init__(self, objects=colors)

    def is_compatible(self, col):
        if col is None:
            return False
        if self._type is None:
            self._type = col[0].__class__.__name__
        if self._type == col[0].__class__.__name__:
            return True
        else:
            return False

    def update(self):
        self.com.pyfem.set_ncolors(self.com.ncolors)
        for Obj in self._list:
            Obj.update()

class manager_operators(manager):
    def __init__(self, colors=[]):
        self._type = "oper"
        manager.__init__(self, colors=colors)

    def update(self):
        manager.update(self)

class manager_fields(manager):
    def __init__(self, colors=[]):
        self._type = "field"
        manager.__init__(self, colors=colors)

    def update(self):
        manager.update(self)

if __name__ == '__main__':
    # --------------------------------------
    # basic class test
    # --------------------------------------
    print(">>> basic class test")
    blue = color([2,22,222])
    blue.append(1)
    blue.append(11)
    blue.append(111)
    print("blue.n = ", blue.n)
    for obj in blue:
        print(obj)
    # --------------------------------------

    # --------------------------------------
    # operator class test
    # --------------------------------------
    print(">>> operator class test")
    from pigasus.gallery.poisson import poisson
    from pigasus.gallery.bilaplacian import bilaplacian
    from caid.cad_geometry import square
    geo = square(n=[3,3], p=[2,2])
    PDE_1 = bilaplacian(geometry=geo)
    PDE_2 = poisson(geometry=geo, V=PDE_1.V)
    PDE_3 = poisson(geometry=geo, V=PDE_1.V)
    PDE_4 = poisson(geometry=geo)

#    print id(PDE_1.V), PDE_1.V.id
#    print id(PDE_2.V), PDE_2.V.id
#    print "---"
#    print PDE_1.operators
#    print PDE_2.operators

    S_1 = PDE_1.D2
    S_2 = PDE_2.stiffness
    S_3 = PDE_3.stiffness
    S_4 = PDE_4.stiffness

    green = color_operator()
    red   = color_operator()

    red.append(S_1)
    green.append(S_2)
    green.append(S_3)
    green.append(S_4)
    print("green.n = ", green.n)
    for obj in green:
        print(id(obj))
    # --------------------------------------

    # --------------------------------------
    # field class test
    # --------------------------------------
    print(">>> field class test")

    F_1 = PDE_1.rhs
    F_2 = PDE_2.rhs
    U_3 = PDE_3.unknown
    U_4 = PDE_4.unknown

    white = color_field()
    white.append(F_1)
    white.append(F_2)
    white.append(U_3)
    white.append(U_4)

    print("white.n = ", white.n)
    for obj in white:
        print(id(obj))
    # --------------------------------------


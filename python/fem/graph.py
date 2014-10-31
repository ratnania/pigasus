# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__all__ = ['graph']

from . import common_obj as _com
from . import constants as _cst
import numpy as _np
from .pigasusObject import *

class graph(pigasusObject):
    def __init__(self, spaces=None):
        pigasusObject.__init__(self)

        self.spaces = spaces

        if spaces is None:
            """
            this means that the current graph was created by bolck
            pigasus will create and compute the graph and then set it
            directly into fortran
            """
            self.spaces = [-1,-1]
            print("graph: Not yet implemented")
            import sys; sys.exit(0)


        # this must be the last thing to do
        self.id = self.com.ngraphs
        self.com.ngraphs += 1
        self.com.graphs.append(self)

        self.setInfoData()

    def setInfoData(self):
        """
        prints informations about the current matrix
        """
        self.infoData['id'] = str(self.id)
        self.infoData['spaces'] = str([V.id for V in self.spaces])



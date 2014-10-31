# -*- coding: UTF-8 -*-
#! /usr/bin/python

__author__="ARA"
__all__ = ['tracelog']
__date__ ="$Mai 09, 2014 10:34:00 PM$"

from . import common_obj as com
class tracelog:
    def __init__(self):
        self._com = com.common_obj()

    def printlog(self, message, condition=True, level=0):
        self._com.pyfem.pyfem_printlog(message, condition, level)

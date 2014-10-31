# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__date__ ="$Jan 12, 2012 3:48:48 PM$"

def _singleton(cls):
    instances = {} # Line 2
    def getinstance():
        if cls not in instances:
            instances[cls] = cls() # Line 5
        return instances[cls]
    return getinstance

@_singleton
class common_obj(object):

    def __init__(self):
        try:
            import core as co
            self.pyfem = co.pyfem
            self.initialize()
        except:
            print("Error while importing pyfem. Pigasus will stop immediatly")
            raise

    def reset(self):
        deleted = True
        self.initialize()

    def initialize(self):
        self.usefiga = False

        self.ngraphs         = 0
        self.nfields         = 0
        self.nmatrices       = 0
        self.noperators      = 0
        self.nspaces         = 0
        self.nmappings       = 0
        self.nmetrics        = 0
        self.ngrids          = 0
        self.nnorms          = 0
        self.nsolvers        = 0

        # ... lists
        self.fields         = []
        self.matrices       = []
        self.operators      = []
        self.spaces         = []
        self.mappings       = []
        self.metrics        = []
        self.grids          = []
        self.norms          = []
        self.graphs         = []
        self.solvers        = []

        # for each grids, we store the id of related fields
        self.list_Grfields_id = [[]]
        self.list_Grdfields_id = [[]]
        # for each grids, we store the id of related matrices
        self.list_Groperators_id = [[]]
        self.list_Grnorms_id = [[]]


def isFloat(Obj):
    return Obj.__class__.__name__=="float"

def isList(Obj):
    return Obj.__class__.__name__=="list"

def isNumpyArray(Obj):
    return Obj.__class__.__name__=="ndarray"

def isField(Obj):
    return Obj.__class__.__name__=="field"

def isMatrix(Obj):
    return Obj.__class__.__name__=="matrix"

def isScipyMatrix(Obj):
    return Obj.__class__.__name__ in ["lil_matrix" , "coo_matrix" , "csr_matrix" , "csc_matrix"]

# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ARA"
__date__ ="$Dec 21, 2011 8:48:23 AM$"

from connectivity import *
from numpy import array

class connectivity_1D(connectivity):
    def __init__(self, geometry, ai_ndof=1):
        connectivity.__init__(self, geometry, ai_ndof )
        self.dim = 1

        self.init_ELT(geometry)

    def init_ELT(self,geometry):
        for li_id in range(0, self.npatch):
            nrb = geometry[li_id]
            list_info_knots = self.list_info_knots(nrb)

            list_Elt_Index = []
            # acces a l 'indice du noeud
            for i in list_info_knots[0][2]:
                list_Elt_Index.append([i])

            self.list_Elt_Index.append(list_Elt_Index)

    def init_data_structure(self, bound_cond):

        for li_id in range(0, self.npatch):

            lpi_Elt_Index = self.list_Elt_Index[li_id]
            li_nen = self.list_nen[li_id]
            li_nel = self.list_nel[li_id]
            lpi_N = self.list_n[li_id]
            lpi_P = self.list_p[li_id]

            _list_elt_st = self.list_elt_st[li_id]
            list_elt_st = [ [i for i in _list_elt_st[0] if lpi_P[0]+1<=i and i<=lpi_N[0]]]
#            print "self.list_elt_st=", self.list_elt_st
#            print "list_elt_st=", list_elt_st
#            print "li_nel=", li_nel

            import numpy as np
            list_IEN  = np.zeros((li_nen, li_nel), dtype=np.int)

            li_baseA = 0
            li_baseB = 0

            for li_dof in range(0, self.ndof):
                li_e = 0
                for li_i in list_elt_st[0]:

                    # A starts from 1
                    li_A = li_i
                    li_A = li_A + li_baseA
                    for li_iloc in range(0, lpi_P[0] + 1):

                        li_B = li_A -(lpi_P[0] - li_iloc)
                        li_Bloc = li_iloc
                        li_Bloc = li_Bloc + li_baseB

                        # as A starts from 1, we must deduce 1
                        list_IEN[li_Bloc, li_e] = ( li_B - 1 )

                    li_e += 1

                li_baseA = li_baseA + li_A
                li_baseB = li_baseB + li_nen

            self.IEN.append(list_IEN)

        # INITIALIIZING THE ID ARRAY
        self.init_ID(bound_cond)

        # INITIALIIZING THE LOCATION MATRIX LM ARRAY
        self.init_LM()

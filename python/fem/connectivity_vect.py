# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="root"
__date__ ="$Mar 28, 2012 11:47:08 AM$"

if __name__ == "__main__":
    print("Hello World");

from .connectivity import *

class connectivity_vect(connectivity):
    def __init__(self, spaces, ai_ndof):
        connectivity.__init__(self, spaces[0].geometry, ai_ndof )

        self.init_ELT(spaces)

    def init_ELT(self,spaces):
#        print "spaces[0].connectivity.list_Elt_Index=", spaces[0].connectivity.list_Elt_Index
        self.list_Elt_Index = spaces[0].connectivity.list_Elt_Index

    def init_ID(self,spaces):
        li_N   = sum ( [len(S.connectivity.ID) for S in spaces] )
        li_npatch = spaces[0].connectivity.npatch

        import numpy as np
        self.ID = np.zeros(li_N, dtype=np.int)

        # initiliazing ID
        li_base_N   = 0
        li_base_ID  = 0
        for S in spaces :
            li_Ni = len(S.connectivity.ID)
#            print "li_Ni=", li_Ni
#            print "li_base_ID=", li_base_ID
#            self.ID[li_base_N : li_base_N+li_Ni] = li_base_ID + S.connectivity.ID[:]
            for i in range(0,li_Ni) :
                if S.connectivity.ID[i] != 0:
                    self.ID[li_base_N+i] = li_base_ID + S.connectivity.ID[i]

            li_base_ID += max (S.connectivity.ID)
            li_base_N  += li_Ni

#        print "self.ID=", self.ID

    def init_data_structure(self,spaces):
        li_N   = sum ( [len(S.connectivity.ID) for S in spaces] )
        li_npatch = spaces[0].connectivity.npatch

        import numpy as np

        # initiliazing IEN
        li_maxnen = 0
        li_maxnel = 0

        li_base_nen = 0
        for li_id in range(0,li_npatch):
            li_nen = 0
            for S in spaces :
                li_nen += S.connectivity.list_nen[li_id]
            if li_nen > li_maxnen :
                li_maxnen = li_nen

            li_nel  = spaces[0].connectivity.list_nel[li_id]
            if li_nel > li_maxnel :
                li_maxnel = li_nel

            lpi_IEN  = np.zeros((li_nen,li_nel), dtype=np.int)

            li_base_nen = 0
            li_base_nnp = 0
            for S in spaces :
                li_neni = S.connectivity.list_nen[li_id]
                lpi_IEN [li_base_nen:li_base_nen+li_neni,:] = S.connectivity.IEN[li_id][:,:]
                lpi_IEN [li_base_nen:li_base_nen+li_neni,:] += li_base_nnp
                li_base_nnp += S.connectivity.list_nnp[li_id]
                li_base_nen += li_neni

            self.IEN.append(lpi_IEN)

#        print "add_connectivity for space_vect"
#        print "li_N, li_npatch, li_maxnen, li_maxnel=", li_N, li_npatch, li_maxnen, li_maxnel


        self.maxnel = li_maxnel
        self.maxnen = li_maxnen

        self.size = sum ( [S.size for S in spaces] )

        # INITIALIIZING THE ID ARRAY
        self.init_ID(spaces)

#        print "self.ID =", self.ID
#        print "self.IEN=", self.IEN

        # INITIALIIZING THE LOCATION MATRIX LM ARRAY
        self.init_LM()

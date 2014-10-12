# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="abdelkaderratnani"
__date__ ="$Jan 20, 2012 7:41:03 PM$"

from boundary_conditions import *
class boundary_conditions_1D(boundary_conditions):
    def all_dirichlet(self, geometry, list_interior_ind=[[]]):
        """
        bc = 10 with the old version of PyIGA
        we must exlude interior boundaries between 2 (or more) patchs,
        this is given in list_interior_ind
        """
        self.list_Dirichlet_ind = []

        # computing the number of patchs
        li_npatch = geometry.npatchs

        for li_id in range(0, li_npatch):
            lo_domain = geometry[li_id]

            list_n = lo_domain.shape

            list_Dirichlet_ind = []

            li_i = 0
            if [li_i] not in list_interior_ind :
                list_Dirichlet_ind.append([li_i])

            li_i = list_n[0] - 1
            if [li_i] not in list_interior_ind :
                list_Dirichlet_ind.append([li_i])

            self.list_Dirichlet_ind.append(list_Dirichlet_ind)

#            print "list_Dirichlet_ind=", list_Dirichlet_ind

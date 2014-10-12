# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="root"
__date__ ="$Dec 22, 2011 10:57:22 AM$"

from boundary_conditions import *
class boundary_conditions_2D(boundary_conditions):
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

            for li_i in range (0, list_n[0]):
                li_j = 0
                if [li_i,li_j] not in list_interior_ind :
                    list_Dirichlet_ind.append([li_i,li_j])

                li_j = list_n[1] - 1
                if [li_i,li_j] not in list_interior_ind :
                    list_Dirichlet_ind.append([li_i,li_j])

            for li_j in range (0, list_n[1]):
                li_i = 0
                if [li_i,li_j] not in list_interior_ind :
                    list_Dirichlet_ind.append([li_i,li_j])

                li_i = list_n[0] - 1
                if [li_i,li_j] not in list_interior_ind :
                    list_Dirichlet_ind.append([li_i,li_j])

            self.list_Dirichlet_ind.append(list_Dirichlet_ind)

            print "list_Dirichlet_ind=", list_Dirichlet_ind

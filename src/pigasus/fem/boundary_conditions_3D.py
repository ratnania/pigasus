# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="ratnani"
__date__ ="$Jan 22, 2012 2:26:39 PM$"

from boundary_conditions import *
class boundary_conditions_3D(boundary_conditions):
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

            for li_k in range (0, list_n[2]):
                for li_j in range (0, list_n[1]):
                    li_i = 0
                    if [li_i,li_j,li_k] not in list_interior_ind :
                        list_Dirichlet_ind.append([li_i,li_j,li_k])

                    li_i = list_n[0] - 1
                    if [li_i,li_j,li_k] not in list_interior_ind :
                        list_Dirichlet_ind.append([li_i,li_j,li_k])

            for li_j in range (0, list_n[1]):
                for li_i in range (0, list_n[0]):
                    li_k = 0
                    if [li_i,li_j,li_k] not in list_interior_ind :
                        list_Dirichlet_ind.append([li_i,li_j,li_k])

                    li_k = list_n[2] - 1
                    if [li_i,li_j,li_k] not in list_interior_ind :
                        list_Dirichlet_ind.append([li_i,li_j,li_k])

            for li_i in range (0, list_n[0]):
                for li_k in range (0, list_n[2]):
                    li_j = 0
                    if [li_i,li_j,li_k] not in list_interior_ind :
                        list_Dirichlet_ind.append([li_i,li_j,li_k])

                    li_j = list_n[1] - 1
                    if [li_i,li_j,li_k] not in list_interior_ind :
                        list_Dirichlet_ind.append([li_i,li_j,li_k])

            self.list_Dirichlet_ind.append(list_Dirichlet_ind)

#            print "list_Dirichlet_ind=", list_Dirichlet_ind

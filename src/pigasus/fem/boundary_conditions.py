# -*- coding: UTF-8 -*-
#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="root"
__date__ ="$Dec 22, 2011 9:42:24 AM$"
__all__ ="boundary_conditions"

def indices_face(nrb, ai_face, shift):
    li_dim = nrb.dim
    if li_dim == 1:
        return _indices_face_1D(nrb, ai_face, shift)
    if li_dim == 2:
        return _indices_face_2D(nrb, ai_face, shift)
    if li_dim == 3:
        return _indices_face_3D(nrb, ai_face, shift)

def _indices_face_1D(nrb, ai_face, shift):
    """
    this routine return the list of indices that are on the face n ai_face
    """
    list_n = nrb.shape

    list_ind = []

    # face n : 1
    if ai_face == 0:
        li_i = 0 + shift
        list_ind.append([li_i])
        return list_ind

    # face n : 2
    if ai_face == 1:
        li_i = list_n[0] - 1 + shift
        list_ind.append([li_i])
        return list_ind

    print "Error _indices_face_1D: you gave a wrong face id. Given face ", ai_face

    import sys
    sys.exit(1)

def _indices_face_2D(nrb, ai_face, shift):
    """
    this routine return the list of indices that are on the face n ai_face
    """
    list_n = nrb.shape

    list_ind = []

    # face n : 1
    if ai_face == 0:
        for li_i in range (0, list_n[0]):
            li_j = 0 + shift
            list_ind.append([li_i,li_j])
        return list_ind

    # face n : 3
    if ai_face == 2:
        for li_i in range (0, list_n[0]):
            li_j = list_n[1] - 1 + shift
            list_ind.append([li_i,li_j])
        return list_ind

    # face n : 2
    if ai_face == 1:
        for li_j in range (0, list_n[1]):
            li_i = 0 + shift
            list_ind.append([li_i,li_j])
        return list_ind

    # face n : 4
    if ai_face == 3:
        for li_j in range (0, list_n[1]):
            li_i = list_n[0] - 1 + shift
            list_ind.append([li_i,li_j])
        return list_ind

    print "Error _indices_face_2D: you gave a wrong face id. Given face ", ai_face
    import sys
    sys.exit(1)

def _indices_face_3D(nrb, ai_face, shift):
    """
    this routine return the list of indices that are on the face n ai_face
    """
    list_n = nrb.shape

    list_ind = []

    # face n : 1
    if ai_face == 0:
        for li_k in range (0, list_n[2]):
            for li_j in range (0, list_n[1]):
                li_i = 0 + shift
                list_ind.append([li_i,li_j,li_k])
        return list_ind

    # face n : 4
    if ai_face == 3:
        for li_k in range (0, list_n[2]):
            for li_j in range (0, list_n[1]):
                li_i = list_n[0] - 1 + shift
                list_ind.append([li_i,li_j,li_k])
        return list_ind

    # face n : 2
    if ai_face == 1:
        for li_j in range (0, list_n[1]):
            for li_i in range (0, list_n[0]):
                li_k = 0 + shift
                list_ind.append([li_i,li_j,li_k])
        return list_ind

    # face n : 5
    if ai_face == 4:
        for li_j in range (0, list_n[1]):
            for li_i in range (0, list_n[0]):
                li_k = list_n[2] - 1 + shift
                list_ind.append([li_i,li_j,li_k])
        return list_ind

    # face n : 3
    if ai_face == 2:
        for li_i in range (0, list_n[0]):
            for li_k in range (0, list_n[2]):
                li_j = 0 + shift
                list_ind.append([li_i,li_j,li_k])
        return list_ind

    # face n : 6
    if ai_face == 5:
        for li_i in range (0, list_n[0]):
            for li_k in range (0, list_n[2]):
                li_j = list_n[1] - 1 + shift
                list_ind.append([li_i,li_j,li_k])
        return list_ind

    print "Error : you gave a wrong face id"
    import sys
    sys.exit(1)

class boundary_conditions(object):
    def __init__(self, geometry, list_Dirichlet_ind = [], list_Neumann_ind = [] \
    , list_Periodic_ind = [], list_duplicated_ind = [], list_duplicata_ind = []):
        self.DirFaces               = [[]]*geometry.npatchs
        self.DuplicatedFaces        = []
        self.DuplicataFaces         = []
        self.list_Dirichlet_ind     = list_Dirichlet_ind
        self.list_Neumann_ind       = list_Neumann_ind
        self.list_Periodic_ind      = list_Periodic_ind
        self.list_duplicated_ind    = list_duplicated_ind
        self.list_duplicata_ind     = list_duplicata_ind
        # initialization with respect to the number of patchs
        li_npatch = geometry.npatchs
        list_empty = [[]] * li_npatch
        self.dirichlet(geometry, list_empty)
        self.duplicate(geometry, faces_base=None, faces=None)

    def dirichlet(self, geometry, faces):
        """
        we must exlude interior boundaries between 2 (or more) patchs,
        this is given in list_interior_ind
        for each patch we must give the list of the correspondant faces
        on which we would like to put the boundary condition
        """
        self.DirFaces = faces

        self.list_Dirichlet_ind = []

        # computing the number of patchs
        li_npatch = geometry.npatchs

        for li_id in range(0, li_npatch):
            nrb = geometry[li_id]
            list_faces = faces[li_id]
            list_Dirichlet_ind = []
            list_shift = [0] * len(list_faces)

            for (F, Sh) in zip(list_faces, list_shift):
                list_indices = indices_face(nrb, F, Sh)
                for Ind in list_indices:
                    list_Dirichlet_ind.append(Ind)

            self.list_Dirichlet_ind.append(list_Dirichlet_ind)

    def duplicate(self, geometry, faces_base, faces, shift_base=None, shift=None):
        """
        this routine is used to duplicate the basis fct
        faces_base and faces, must be of the form : list of [patch_id, face_id, shift]
        """

        if faces_base is None:
            li_npatch = geometry.npatchs
            self.list_duplicated_ind = [[]] * li_npatch
            self.list_duplicata_ind  = [[]] * li_npatch
            return
        else:
            self.DuplicatedFaces = faces_base
            self.DuplicataFaces = faces

        # verifier que l'indice de face_base est < a celui de face
        for (pf_base, pf) in zip(faces_base, faces):
            pbase_id = pf_base[0] ; fbase_id = pf_base[1]
            p_id = pf[0] ; f_id = pf[1]
            if ( pbase_id == p_id ) and (fbase_id > f_id) :
                print "Error indices for faces: fbase_id must be smaller than or equal to f_id"
                import sys; sys.exit(0)

        self.list_duplicated_ind = []
        self.list_duplicata_ind  = []
        for (pf_base, pf) in zip(faces_base, faces):
#            print "pf_base =", pf_base, "       pf =", pf
            pbase_id = pf_base[0] ; fbase_id = pf_base[1]
            p_id = pf[0] ; f_id = pf[1]
            # if shift has been given
            if len(pf) == 2:
                shbase = 0
                sh = 0
            else:
                shbase = pf_base[2]
                sh = pf[2]

            # getting domains
            nrb_base = geometry[pbase_id]
            nrb       = geometry[p_id]

            list_ind_base = indices_face(nrb_base, fbase_id,shbase)
            list_indices  = indices_face(nrb, f_id,sh)
            for (Ind_base, Ind) in zip(list_ind_base, list_indices):
#                print "Ind_base = ", Ind_base, "     Ind = ", Ind
                list_Indbase = [pbase_id]
                for I in Ind_base:
                    list_Indbase.append(I)
                list_Ind = [p_id]
                for I in Ind:
                    list_Ind.append(I)
                self.list_duplicated_ind.append(list_Indbase)
                self.list_duplicata_ind.append(list_Ind)
#        print "self.list_duplicated_ind=", self.list_duplicated_ind
#        print "self.list_duplicata_ind=", self.list_duplicata_ind

#        # computing the number of patchs
#        li_npatch = len(geometry)
#
#        for li_id in range(0, li_npatch):
#            nrb = geometry[li_id]
#            list_faces_base = faces_base[li_id]
#            list_faces      = faces[li_id]
#            if shift_base is None:
#                list_shift_base = [0] * len(list_faces_base)
#            else:
#                list_shift_base = shift_base[li_id]
#            if shift is None:
#                list_shift = [0] * len(list_faces)
#            else:
#                list_shift      = shift[li_id]
#            list_duplicata_ind  = []
#            list_duplicated_ind = []
#
#            for (F_base, F, Sh_base, Sh) in zip(list_faces_base,list_faces,list_shift_base,list_shift):
#                list_ind_base = nrb.indices_face(F_base,Sh_base)
#                list_indices  = nrb.indices_face(F,Sh)
#                for (Ind_base, Ind) in zip(list_ind_base, list_indices):
#                    list_duplicated_ind.append(Ind_base)
#                    list_duplicata_ind.append(Ind)
#
#            self.list_duplicata_ind.append(list_duplicata_ind)
#            self.list_duplicated_ind.append(list_duplicated_ind)

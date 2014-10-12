# -*- coding: UTF-8 -*-
#! /usr/bin/python

__author__ = "ratnani"
__date__ = "$Dec 20, 2011 3:56:49 PM$"
__all__ = ['connectivity']

import numpy as np
array = np.array
class connectivity:
    def __init__(self, geometry,ai_ndof=1):

        # computing the number of patchs
        self.npatch = geometry.npatchs

        # for each patch, we copy n,p, nel, nen
        self.list_n = []
        self.list_p = []
        self.list_nel = []
        self.list_nen = []
        self.list_nnp = []
        self.list_elt_st = [] # the left-index of the element, per direction
        for li_id in range(0, self.npatch):
            nrb = geometry[li_id]
            li_dim = nrb.dim

            self.list_n.append(nrb.shape)
            self.list_p.append(nrb.degree)

            list_info_knots = self.list_info_knots(nrb)

            # calcul des nen
            li_nen = 1
            for li_d in range(0, li_dim):
                li_nen *= nrb.degree[li_d] + 1
            self.list_nen.append(li_nen)

            # indexation des elements
            li_nel = 1
            for li_d in range(0, li_dim):
                # we test if the knot vector is periodic
                li_n = nrb.shape[li_d]
                li_p = nrb.degree[li_d]
                if len(list_info_knots[li_d][0]) == li_n + li_p + 1:
                    li_nel *= li_n - li_p
                else :
                    li_nel *= list_info_knots[li_d][0].__len__() - 1
            self.list_nel.append(li_nel)

            # calcul des nnp
            li_nnp = 1
            for li_d in range(0, li_dim):
                # we test if the knot vector is periodic
                li_n = nrb.shape[li_d]
                li_p = nrb.degree[li_d]
                if len(list_info_knots[li_d][0]) == li_n + li_p + 1:
                    li_nnp *= li_n
                else:
                    li_nnp *= nrb.shape[li_d]
            self.list_nnp.append(li_nnp)

            # compute list_elt_st
            list_elt_st = []
            for li_d in range(0, li_dim):
                list_elt_st.append(list_info_knots[li_d][2][0:-1])
            self.list_elt_st.append(list_elt_st)

        self.dim = geometry.dim
        self.ndof = ai_ndof

        self.nnp = sum (self.list_nnp)
        self.size = 0
        self.maxnel = 0
        self.maxnen = 0

        self.baseID = [] # this is the length of ID for each patch

        self.ID  = []
        self.IEN = []
        self.LM  = []
        self.ID_loc  = [[]]*self.npatch

        self.list_Elt_Index = []

    def init_ID(self, bound_cond):
        list_n = self.list_n

        DirFaces = bound_cond.DirFaces
        DuplicatedFaces = bound_cond.DuplicatedFaces
        DuplicataFaces = bound_cond.DuplicataFaces

#        print " DirFaces ", DirFaces
#        print " DuplicataFaces ", DuplicataFaces
#        print "DuplicatedFaces  ", DuplicatedFaces

        from idutils import computeLocalID, computeGlobalID
        list_id = computeLocalID(list_n, DirFaces, DuplicatedFaces, DuplicataFaces)
        ID = computeGlobalID(list_id)

        self.ID_loc = list_id
        self.ID = array(ID)

        self.size = max(self.ID)

    def init_LM(self):
#        print "self.baseID=", self.baseID
        import numpy as np
        for li_id in range(0, self.npatch):
            size = self.ID_loc[li_id].size
            ID_loc = list(self.ID_loc[li_id].transpose().reshape(size))
#            print "li_id = ", li_id
#            print "self.ID=", self.ID
#            print "self.IEN[li_id]=", self.IEN[li_id]
#            list_LM = li_maxLM + np.asarray(self.ID[self.IEN[li_id][:,:]])
#            print "self.baseID = ", self.baseID
#            li_baseID = self.baseID[li_id]
#            print "li_baseID = ", li_baseID
            list_LM = []
            lpr_IEN = np.asarray(self.IEN[li_id]).transpose()
            for P in lpr_IEN:
                Q = []
#                print "P = ", P
                for t in P:
#                    print "t, ID_loc[t]=", t, ID_loc[t]
                    Q.append(ID_loc[t])
                list_LM.append(Q)

#            lpi_shape = np.asarray(self.IEN).shape
#            list_maxLM = np.ones(lpi_shape, dtype=np.int)
#            print "li_maxLM = ", li_maxLM
#            list_maxLM *= li_maxLM
#            list_IEN = ( list_maxLM + np.asarray(self.IEN[li_id]).tolist())
#            print "list_maxLM =", list_maxLM
#            list_LM = self.ID[list_IEN]

#            print "list_LM=", list_LM
            self.LM.append(np.asarray(list_LM).transpose())
#        print "self.LM=", self.LM

    def printinfo(self):
        print "*******************************"
        print " global informations "
        print "*******************************"
        print " number of patchs :", self.npatch
        print " nnp :", self.nnp
        print "*******************************"
        for li_id in range(0, self.npatch):
            print "*******************************"
            print " Current Patch-id =", li_id
            li_nel = self.list_nel[li_id]
            li_nen = self.list_nen[li_id]
            li_nnp = self.list_nnp[li_id]
            print "-- nel =", li_nel
            print "-- nen =", li_nen
            print "-- nnp =", li_nnp
            print "*******************************"

            print "======"
            print " IEN "
            print "======"
            for li_e in range(0, li_nel):
                print "     elt =", li_e
                print "          ", self.IEN[li_id][:, li_e]
            print "======"

            print "======"
            print " LM "
            print "======"
            for li_e in range(0, li_nel):
                print "     elt =", li_e
                print "          ", self.LM[li_id][:, li_e]
            print "======"

            print "*******************************"

        print "==============================="
        print " ID "
        print "==============================="
        print " ", self.ID[:]
        print "==============================="

    def save(self, etiq="", fmt='zip', name=None):
        """
        name needed for zip format
        """
        # ...
        def exportTXT(etiq):
            np.savetxt(etiq+"ID.txt"\
                       ,np.asarray(self.ID)\
                      , fmt='%d')
            for li_id in range(0, self.npatch):
                li_nel = self.list_nel[li_id]
                li_nen = self.list_nen[li_id]
                li_nnp = self.list_nnp[li_id]
                np.savetxt(etiq+"IEN_"+str(li_id)+".txt"\
                           ,np.asarray(self.IEN[li_id][:, :])\
                          , fmt='%d')
                np.savetxt(etiq+"LM_"+str(li_id)+".txt"\
                           ,np.asarray(self.LM[li_id][:, :])\
                          , fmt='%d')
        # ...

        # ...
        if fmt == "txt":
            exportTXT(etiq)
        # ...

        # ...
        if (fmt == "zip") and (name is not None):
            import os
            from contextlib import closing
            from zipfile import ZipFile, ZIP_DEFLATED

            # ...
            def zipdir(basedir, archivename):
                assert os.path.isdir(basedir)
                with closing(ZipFile(archivename, "w", ZIP_DEFLATED)) as z:
                    for root, dirs, files in os.walk(basedir):
                        #NOTE: ignore empty directories
                        for fn in files:
                            absfn = os.path.join(root, fn)
                            zfn = absfn[len(basedir)+len(os.sep):] #XXX: relative path
                            z.write(absfn, zfn)
            # ...

            os.system("mkdir -p " + name)
            etiq = name+"/"
            exportTXT(etiq=etiq)
            basedir = name
            archivename = name+".zip"
            zipdir(basedir, archivename)
            os.system("rm -R " + name)
        # ...



    def list_info_knots(self, nrb):
        # contient la multiplicite de chaque noeud,
        # et l'indice du dernier noeud duplique, selon chaque direction
        from itertools import groupby
        li_dim = nrb.dim
        list_info = []
        for li_d in range(0, li_dim):
            lpr_knots = nrb.knots[li_d]
            u = [k for k, g in groupby(lpr_knots)]
            m = [len(list(g)) for k, g in groupby(lpr_knots)]
            index = [sum(m[0:i+1]) for i in range(0,m.__len__())]
            list_info.append((u,m,index))
        return list_info





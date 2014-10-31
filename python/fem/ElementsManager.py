# -*- coding: UTF-8 -*-
#! /usr/bin/python

__author__ = ['ARA']
__all__ = ['ElementsManager']

import numpy as np

# ...
def get_elements1D(nrb):
    dim = nrb.dim

    list_nr = []
    for axis in range(0,dim):
        t = nrb.knots[axis]
        p = nrb.degree[axis]
        ut = np.unique(t) # unique knots
        nut = len(ut)
        nt = len(t)
        nr = nt - nut - p
        list_nr.append(nr)

    nx = nrb.shape[0] ; px = nrb.degree[0] ; tx = nrb.knots[0]
    eltInd = np.zeros(nx+px, dtype = np.int)
    b = 0
    for i in range(0, nx + px):
        meas = ( tx[i+1] - tx[i] )
        if np.abs(meas) < 1.e-12:
            eltInd[i] = -1
        else:
            eltInd[i] = b
            b += 1

    # find for each spline the elements where it lives
    list_supports = []
    for i in range(0, nx):
        list_i = list(range(i, i + px + 1))

        list_elt = []
        for _i in list_i:
            e = eltInd[_i]
            if e > -1:
                list_elt.append(e)

        list_supports.append(list_elt)

    return list_supports
# ...

# ...
def get_elements2D(nrb):
    dim = nrb.dim

    list_nr = []
    for axis in range(0,dim):
        t = nrb.knots[axis]
        p = nrb.degree[axis]
        ut = np.unique(t) # unique knots
        nut = len(ut)
        nt = len(t)
        nr = nt - nut - p
        list_nr.append(nr)

    nx = nrb.shape[0] ; px = nrb.degree[0] ; tx = nrb.knots[0]
    ny = nrb.shape[1] ; py = nrb.degree[1] ; ty = nrb.knots[1]
    eltInd = np.zeros((nx+px,ny+py), dtype = np.int)
    b = 0
    for j in range(0, ny + py):
        for i in range(0, nx + px):
            meas =    ( ty[j+1] - ty[j] ) \
                    * ( tx[i+1] - tx[i] )
            if np.abs(meas) < 1.e-12:
                eltInd[i,j] = -1
            else:
                eltInd[i,j] = b
                b += 1

    # find for each spline the elements where it lives
    list_supports = []
    for j in range(0, ny):
        for i in range(0, nx):
            list_i = list(range(i, i + px + 1))
            list_j = list(range(j, j + py + 1))

            list_elt = []
            for _j in list_j:
                for _i in list_i:
                    e = eltInd[_i,_j]
                    if e > -1:
                        list_elt.append(e)

            list_supports.append(list_elt)

    return list_supports
# ...

# ...
def get_elements3D(nrb):
    dim = nrb.dim

    list_nr = []
    for axis in range(0,dim):
        t = nrb.knots[axis]
        p = nrb.degree[axis]
        ut = np.unique(t) # unique knots
        nut = len(ut)
        nt = len(t)
        nr = nt - nut - p
        list_nr.append(nr)

    nx = nrb.shape[0] ; px = nrb.degree[0] ; tx = nrb.knots[0]
    ny = nrb.shape[1] ; py = nrb.degree[1] ; ty = nrb.knots[1]
    nz = nrb.shape[2] ; pz = nrb.degree[2] ; tz = nrb.knots[2]
    eltInd = np.zeros((nx+px,ny+py,nz+pz), dtype = np.int)
    b = 0
    for k in range(0, nz + pz):
        for j in range(0, ny + py):
            for i in range(0, nx + px):
                meas =    ( tz[k+1] - tz[k] ) \
                        * ( ty[j+1] - ty[j] ) \
                        * ( tx[i+1] - tx[i] )
                if np.abs(meas) < 1.e-12:
                    eltInd[i,j,k] = -1
                else:
                    eltInd[i,j,k] = b
                    b += 1

    # find for each spline the elements where it lives
    list_supports = []
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                list_i = list(range(i, i + px + 1))
                list_j = list(range(j, j + py + 1))
                list_k = list(range(k, k + pz + 1))

                list_elt = []
                for _k in list_k:
                    for _j in list_j:
                        for _i in list_i:
                            e = eltInd[_i,_j,_k]
                            if e > -1:
                                list_elt.append(e)

                list_supports.append(list_elt)

    return list_supports
# ...

# ...
def get_elements(geometry):
    dim = geometry.dim
    list_elts = []
    for nrb in geometry:
        if dim == 1:
            list_elts.append(get_elements1D(nrb))
        if dim == 2:
            list_elts.append(get_elements2D(nrb))
        if dim == 3:
            list_elts.append(get_elements3D(nrb))
    return list_elts
# ...

# ...
class ElementsManager(object):
    """
    this class is used for parallelization
    it is assumed that the partitioning process has been done with an external
    tool (scotch, murge)
    By calling SetLocalData, it will compute the local patch ids and the
    elements on which we need to call the assembly function.
    After Initialization, the user needs to call GetLocalPatchs and
    GetLocalElements
    """
    class support(object):
        def __init__(self, patch_id, list_elts):
            self.id = patch_id
            self.elts = list_elts

    def __init__(self, geometry, nodelist=None, LocalID=None):
        self.geometry = geometry
        self.supports = get_elements(geometry)
        self.nodelist = nodelist
        self.LocalID  = LocalID
        self.patchs = []
        self.elts   = []
        self.IDsupports = []

    def setLocalData(self, nodelist):
        """
        sets local patchs ids and corresponding elements
        """
        nptachs = self.geometry.npatchs
        for id in range(0, npatchs):
            nrb = self.geometry[id]
            supports = self.supports[id]
        pass

    def setNodeList(self, nodelist):
        self.nodelist = nodelist

    def GetLocalPatchs(self):
        """
        returns local patchs ids
        """
        list_patchs = []
        for I in self.nodelist:
            supports = self.IDsupports[I-1]
            list_patchs += [i for [i,elts] in supports]

        list_patchs = np.unique(np.asarray(list_patchs))
        return list(list_patchs)

    def GetLocalElements(self, id):
        """
        returns local elements ids for the giving patch
        """
        list_elts = []
        for I in self.nodelist:
            supports = self.IDsupports[I-1]
            lelts = [elts for [i,elts] in supports if i==id]
#            print "lets = ", lelts
            for ets in lelts:
                for e in elts:
                    list_elts += [e]
        list_elts = np.unique(np.asarray(list_elts))
        list_elts += 1 # must add 1, because the indexation is 0 based
        return list(list_elts)

    def InitializeIDsupports(self, LocalID, ID):
        dim = self.geometry.dim
        def _tolist(ids):
            if dim == 1:
                return ids.transpose()
            if dim == 2:
                return ids.transpose().reshape(ids.size)
            if dim == 3:
                print("Not yet Implemented")

        self.LocalID = LocalID
        self.ID = ID
        nID = np.max(np.array(ID))
        self.IDsupports = [[] for i in range(0,nID)]
#        print "====="
#        print len(self.IDsupports)
#        print "====="
        for ipatch,list_id in enumerate(LocalID):
            ids = _tolist(list_id)
#            print "-------------"
#            print 'ids : ', ids
            Localsupports = self.supports[ipatch]
            nb = len(ids)
            for b in range(0, nb):
                I = ids[b]
                elts = Localsupports[b]
                if I > 0:
#                    print 'I = ', I
#                    print 'elts = ', elts
#                    print self.IDsupports[I-1]
                    self.IDsupports[I-1].append([ipatch, elts])

# ...

# ---------------------------------------------------------
if __name__ == '__main__':
    def test1D():
        print("Test 1D:")
        from caid.cad_geometry import line as domain
        n = [10]
        p = [2]
        geo = domain(n=n,p=p)
        from time import time
        t_start = time()
        list_supports = get_elements(geo)
        t_end = time()
        print("Elapsed time : ", t_end - t_start)
#        print "Supports ", list_supports

    def test2D():
        print("Test 2D:")
        from caid.cad_geometry import square as domain
        n = [10,10]
        p = [2,2]
        geo = domain(n=n,p=p)
        from time import time
        t_start = time()
        list_supports = get_elements(geo)
        t_end = time()
        print("Elapsed time : ", t_end - t_start)
#        print "Supports ", list_supports

    def test3D():
        print("Test 3D:")
        from caid.cad_geometry import trilinear as domain
        n = [10,10,10]
        p = [2,2,2]
        geo = domain(n=n,p=p)
        from time import time
        t_start = time()
        list_supports = get_elements(geo)
        t_end = time()
        print("Elapsed time : ", t_end - t_start)
#        print "Supports ", list_supports

    def test2D_EM():
        print("Test 2D-EM:")
#        PRINT = True
        PRINT = False
        from caid.cad_geometry import square as domain
        n = [50,50]
        p = [3,3]
        geo = domain(n=n,p=p)
        nrb = geo[0]
        # creation of a multi-patch domaine with 4 patchs
        npatchs = 4
        for i in range(1,npatchs):
            geo.append(nrb)
        # with the following boundary Conditions
        DirFaces = [[1,2],[2,3],[0,3],[0,1]]
        DuplicatedFaces = [[0,3],[1,0],[2,1],[0,0]]
        DuplicataFaces  = [[1,1],[2,2],[3,3],[3,2]]

        list_n = [nrb.shape for nrb in geo]
        from .idutils import computeLocalID, computeGlobalID, print_id
        list_id = computeLocalID(list_n, DirFaces, DuplicatedFaces, DuplicataFaces)
        ID = computeGlobalID(list_id)

        from time import time
        t_start = time()
        EM = ElementsManager(geo)
        list_supports = EM.supports
        EM.InitializeIDsupports(list_id, ID)
        t_end = time()
        print("Elapsed time : ", t_end - t_start)

        t_start = time()
        nodelist = list(range(0, (n[0]+2*p[0]+1)/2))
        EM.setNodeList(nodelist)
        list_patchs = EM.GetLocalPatchs()
#        print "listpatchs = ", list_patchs
        for ipatch in list_patchs:
            list_elts = EM.GetLocalElements(ipatch)
#            print "patch ", i, " elements ", list_elts

        t_end = time()
        print("Elapsed time : ", t_end - t_start)

        if PRINT :
            print("--------------  FINAL ------------------")
            for i,id in enumerate(list_id):
                print("...... patch id : ", i, " ......")
                print_id(id)
                print("...... supports ................")
                print(list_supports[i])
            print("--------------  ID ------------------")
            print(ID)
            print("--------------  I DATA --------------")
            for I,data in enumerate(EM.IDsupports):
#                print "I ", I, " data ", data
                print("<<< I ", I+1, " >>>")
                for [ipatch, elts] in data:
                    print("patch ", ipatch, " elements ", elts)

    if True:
#    if False:
#        test1D()
#        test2D()
#        test3D()
        test2D_EM()

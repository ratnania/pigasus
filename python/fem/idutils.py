# -*- coding: UTF-8 -*-

__author__ = 'ARA'
__all__ = ['computeLocalID', 'computeGlobalID']

import numpy as np


# ...
def initLocalID(faces, n, base):
    dim = len(n)
    if dim == 1:
        return initLocalID_1D(faces, n, base)
    if dim == 2:
        return initLocalID_2D(faces, n, base)

def initLocalID_1D(faces, n, base):
    id =- np.ones(n, dtype=np.int)

    dim = len(n)

    ib = 0 ; ie = n[0]-1

    for f in faces:

        if f==0:
            ib += 1

        if f==1:
            ie -= 1

    for i in range(ib,ie+1):
        A = i - ib
        id[i] = A + base

    id += 1
    base += ie - ib + 1

    return id, base

def initLocalID_2D(faces, n, base):
    id =- np.ones(n, dtype=np.int)

    dim = len(n)

    ib = 0 ; ie = n[0]-1
    jb = 0 ; je = n[1]-1

    for f in faces:

        if f==0:
            jb += 1

        if f==2:
            je -= 1

        if f==1:
            ib += 1

        if f==3:
            ie -= 1

    ne = ie  - ib + 1

#    print "ib,ie = ", ib,ie
#    print "jb,je = ", jb,je

    for j in range(jb,je+1):
        for i in range(ib,ie+1):
            A = ( j - jb ) * ne + i - ib
            id[i,j] = A + base

    id += 1
    base += ( je - jb ) * ne + ie - ib + 1

    return id, base
# ...
def print_id(id):
    dim = len(id.shape)
    if dim ==1:
        print_id_1D(id)
    if dim ==2:
        print_id_2D(id)

def print_id_1D(id):
    id_ = np.zeros_like(id)
    n, = id.shape
    for i in range(0,n):
        id_[i] = id[i]
    print id_.transpose()

def print_id_2D(id):
    id_ = np.zeros_like(id)
    n,m = id.shape
    for j in range(0,m):
        for i in range(0,n):
            id_[i,j] = id[i,-j-1]
    print id_.transpose()
    print id.transpose().reshape(id.size)
# ...

# ...
def isDuplicata(patch_id, face, DuplicataPatchs):
    for data in DuplicataPatchs:
        if (data[0]==patch_id) and (data[1]==face):
            return True
    return False
# ...

# ...
def get_ij_1D(n, f):
    if f==0:
        list_i = [0]

    if f==1:
        list_i = [n[0]-1]

    return list_i

def get_ij_2D(n, f):
    if f==0:
        list_i = range(0, n[0])
        list_j = [0] * n[0]

    if f==1:
        list_i = [0] * n[1]
        list_j = range(0, n[1])

    if f==2:
        list_i = range(0, n[0])
        list_j = [n[1] - 1] * n[0]

    if f==3:
        list_i = [n[0] - 1] * n[1]
        list_j = range(0, n[1])

    return list_i, list_j

# ...
def updateDuplicated_1D(n_m, n_s, list_id, p_m, f_m, p_s, f_s):
    """
    p_m : master patch
    f_m : master face
    p_s : slave patch
    f_s : slave face
    """
    list_i_m = get_ij_1D(n_m, f_m)
    list_i_s = get_ij_1D(n_s, f_s)

    for (i_m,i_s) in zip(list_i_m, list_i_s):
        id_s = list_id[p_s]
        id_m = list_id[p_m]

        id_s[i_s] = id_m[i_m]

    return list_id

def updateDuplicated_2D(n_m, n_s, list_id, p_m, f_m, p_s, f_s):
    """
    p_m : master patch
    f_m : master face
    p_s : slave patch
    f_s : slave face
    """
    list_i_m, list_j_m = get_ij_2D(n_m, f_m)
    list_i_s, list_j_s = get_ij_2D(n_s, f_s)

    for (i_m,j_m,i_s,j_s) in zip(list_i_m, list_j_m, list_i_s, list_j_s):
        id_s = list_id[p_s]
        id_m = list_id[p_m]

        id_s[i_s,j_s] = id_m[i_m,j_m]
#        print ">> ",i_m, j_m, id_m[i_m,j_m]

#    print id_s

    return list_id

def updateDuplicated(n_m, n_s, list_id, p_m, f_m, p_s, f_s):
    dim = len(n_m)
    if dim == 1:
        return updateDuplicated_1D(n_m, n_s, list_id, p_m, f_m, p_s, f_s)
    if dim == 2:
        return updateDuplicated_2D(n_m, n_s, list_id, p_m, f_m, p_s, f_s)
# ...

# ...
def computeLocalID(list_n, DirFaces, DuplicatedFaces, DuplicataFaces):
    dim       = len(list_n[0])
    npatchs   = len(list_n)
    AllFaces  = range(0,2 * dim)
    AllPatchs = range(0,npatchs)

    BasePatchs = [0]
    DuplicatedPatchs = list(np.unique(np.array([data[0] for data in DuplicatedFaces])))
    DuplicataPatchs = list(np.unique(np.array([data[0] for data in DuplicataFaces])))

    base = 0
    list_id = []
    for i in range(0, npatchs):
        list_id.append([])
    for patch_id,faces in enumerate(DirFaces):
        _faces = [f for f in faces]
        # ... mettre a jour faces, en rajoutant les faces dupliquees
        if patch_id in DuplicataPatchs:
            list_faces = [f for f in AllFaces if f not in faces]
            for f in list_faces:
                if isDuplicata(patch_id, f, DuplicataFaces):
                    _faces.append(f)

        id, base = initLocalID(_faces, list_n[patch_id], base)
        list_id[patch_id] = id

#    print "-------------- INIT  ------------------"
#    for i,id in enumerate(list_id):
#        print "...... patch id : ", i, " ......"
#        print_id(id)

    for data_m, data_s in zip(DuplicatedFaces, DuplicataFaces):
        p_m = data_m[0]   ; f_m = data_m[1]
        p_s = data_s[0]   ; f_s = data_s[1]
        n_m = list_n[p_m] ; n_s = list_n[p_s]
        list_id = updateDuplicated(n_m, n_s, list_id, p_m, f_m, p_s, f_s)

    return list_id
# ...

# ...
def computeGlobalID(list_id):
    ID = []
    for id in list_id:
        ID += list(id.transpose().reshape(id.size))
    return ID
# ...

if __name__ == '__main__':
    if False:
        from time import time

        t_start = time()

        PRINT = True
    #    PRINT = False

    ##    list_n = [[4]]*3
    #    list_n = [[1024]]*3
    #    DirFaces = [[0],[],[1]]
    #    DuplicatedFaces = [[0,1],[1,1]]
    #    DuplicataFaces  = [[1,0],[2,0]]

        list_n = [[3,3]]*4
    #    list_n = [[1024,1024]]*4

        DirFaces = [[1,2],[2,3],[0,3],[0,1]]
        DuplicatedFaces = [[0,3],[1,0],[2,1],[0,0]]
        DuplicataFaces  = [[1,1],[2,2],[3,3],[3,2]]

        list_id = computeLocalID(list_n, DirFaces, DuplicatedFaces, DuplicataFaces)
        ID = computeGlobalID(list_id)

        if PRINT :
            print "--------------  FINAL ------------------"
            for i,id in enumerate(list_id):
                print "...... patch id : ", i, " ......"
                print_id(id)
            print "--------------  ID ------------------"
            print(ID)

        t_end = time()
        print "Elapsed time ", t_end - t_start

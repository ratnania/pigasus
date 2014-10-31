# -*- coding: UTF-8 -*-
import numpy as np
from scipy.sparse import eye
from pigasus.fem import common_obj as _com

__author__ = 'ARA'
__all__ = ['BlockVector', 'BlockMatrix']

# ---------------------------------------
def create_block_line(list_A):
    """
     it is supposed that all matrices in list_A have the same number of lines
    """
#    import time
#    t_start = time.time()

    nR = list_A[0].shape[0]
    nC = sum([A.shape[1] for A in list_A])
#    print nR, nC
    nnz = sum([A.nnz for A in list_A])

    ia = np.zeros(nR+1, dtype=int)
    ja = np.zeros(nnz, dtype=int)
    a  = np.zeros(nnz)

    list_nC = np.zeros(len(list_A), dtype=int)
    for i in range(1, len(list_A)):
        list_nC[i] = sum([A.shape[1] for A in list_A[:i]])

    index = 0
    for row in range(0, nR):
        ia[row] = index
        for (A,n) in zip(list_A, list_nC):
            nel = A.indptr[row+1] - A.indptr[row]
            loc_ja = A.indices[A.indptr[row] : A.indptr[row+1]] + n
            ja[index: index + nel] = loc_ja
            a [index: index + nel] = A.data[A.indptr[row] : A.indptr[row+1]]
            index += nel

    ia [nR] = nnz

    from scipy.sparse import csr_matrix
    M = csr_matrix( (a,ja,ia), shape=(nR,nC) )

#    t_end = time.time()
#    print "Elapsed time - Block Line Assembling : ", t_end - t_start

    return M
# ---------------------------------------

# ---------------------------------------
def create_block_column(list_A):
#    import time
#    t_start = time.time()
    list_tA = [A.transpose().tocsr() for A in list_A]
    M = create_block_line(list_tA)
    Mt = M.transpose().tocsr()
#    t_end = time.time()
#    print "Elapsed time - Block Column Assembling : ", t_end - t_start
    return Mt
# ---------------------------------------

# ---------------------------------------
def assembly_block_matrices(matrices):
    list_M = []
    for list_A in matrices:
        list_M.append(create_block_line(list_A))
    M = create_block_column(list_M)
    return M
# ---------------------------------------

# ---------------------------------------

# ---------------------------------------


# ---------------------------------------
class BlockMatrix(object):
    def __new__(typ, *args, **kwargs):
        obj = object.__new__(typ)
        obj._id = None
        return obj

    def __init__(self, matrices):
        self.M = None
        self._id = None

        n = len(matrices)
        m = len(matrices[0])
        self._shape = [n,m]

        # .. find None matrices and replace them by a zero matrix
        list_matrices = []
        for i in range(0,n):
            line = []
            for j in range(0,m):
                M = matrices[i][j]
                if M is None:
                    M_shape = [0,0]
                    # ... find a non-None matrix in the same block line as M
                    for i1 in range(0,n):
                        M1 = matrices[i1][j]
                        if M1 is not None:
                            M_shape[0] = M1.shape[0]
                    # ... find a non-None matrix in the same block column as M
                    for j1 in range(0,m):
                        M1 = matrices[i][j1]
                        if M1 is not None:
                            M_shape[1] = M1.shape[1]

                    Z = 0. * eye(M_shape[0], M_shape[1]) ; Z = Z.tocsr() # zero matrix
                    line.append(Z)
                else:
                    line.append(M)

            list_matrices.append(line)
        self._matrices = list_matrices
        # ...

    @property
    def id(self):
        return self._id

    @property
    def matrices(self):
        return self._matrices

    @property
    def shape(self):
        return self._shape

    def assembly(self):
        list_M = []
        for list_A in self.matrices:
            list_M.append(create_block_line(list_A))
        self.M = create_block_column(list_M)

    def todense(self):
        return self.M.todense()

    def get(self):
        return self.M

    def dot(self, other):
        """
        scalar dot between a matrix and a field
        """
        return self.__mul__(other)

    def __mul__(self, other):
        # ... case of numpy array
        if _com.isNumpyArray(other):
            M = self.get()
            Y = M.dot(other)
            return Y
        # ...
        # ... case of pigasus field
        if _com.isField(other):
            M = self.get()
            Y = M.dot(other.get())
            from pigasus.fem.field import field
            F = field.__new__(field)
            F.set(Y)
            return F
        # ...
        # ... case of list of numpy arrays or pigasus fields
        if _com.isList(other):
            L = []
            for F in other:
                if _com.isNumpyArray(F):
                    L += list(F)
                if _com.isField(F):
                    L += list(F.get())
            X = np.asarray(L)
            M = self.get()
            Y = M.dot(X)
            from pigasus.fem.field import field
            ind_b = 0
            list_F = []
            for U in other:
                F = None
                if _com.isNumpyArray(U):
                    U_size = U.shape[0]
                    F = Y[ind_b:ind_b+U_size]
                    ind_b += U_size
                if _com.isField(U):
                    F = field.__new__(field)
                    F.set(Y[ind_b:ind_b+U.size])
                    ind_b += U.size
                list_F.append(F)
            return list_F
        # ...

#        # ... case of a scalar number
#        if _com.isFloat(other):
#            M = BlockMatrix.__new__(BlockMatrix)
#            M.set(other * self.get ())
#            return M
#        # ...
#        if _com.isScipyMatrix(other):
#            M = BlockMatrix.__new__(BlockMatrix)
#            M.set(other * self.get ())
#            return M
        # ...

        print("Not yet implemented in __mul__")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pos__(self):
        self.__mul__(1.0)

    def __neg__(self):
        self.__mul__(-1.0)

    def __imul__(self, other):
        if _com.isFloat(other):
            self.com.pyfem.matrix_mult_scal (self.id, other, self.id)
        if _com.isMatrix(other):
            self.com.pyfem.matrix_mult_matrix (self.id,other.id,self.id)
        if _com.isScipyMatrix(other):
            self.set(other * self.get ())
        return self
# ---------------------------------------

# ---------------------------------------
class BlockVector(object):
    def __init__(self, vectors):
        self.vectors = vectors
        self.rhs = None

    def assembly(self):
        rhs = []
        for r in self.vectors:
           rhs += list(r)
        self.rhs = np.asarray(rhs)

    def get(self):
        return self.rhs
# ---------------------------------------



if __name__ == '__main__':
    # ---------------------------------------
    def assembly_dense(matrices):
        nR = 0 ; list_nR = [0]
        for list_M in matrices:
            nR += list_M[0].shape[0]
            list_nR.append(nR)
        nC = 0 ; list_nC = [0]
        for M in matrices[0]:
            nC += M.shape[1]
            list_nC.append(nC)
    #    print nR,nC
    #    print list_nR
    #    print list_nC

        Mdense = np.zeros((nR, nC))

        m = len(matrices)
        n = len(matrices[0])

        for j in range(0,n):
            for i in range(0,m):
                Mij = matrices[i][j].todense()
                i_ = list_nR[i] ; j_ = list_nC[j]
                nR_ = Mij.shape[0]; nC_ = Mij.shape[1]
    #            print "nR_, nC_ = ", nR_, nC_
    #            print "i,j, i_,j_ =", i,j,i_,j_
                for J in range(0,nC_):
                    for I in range(0,nR_):
                        Mdense[i_+I, j_+J] = Mij[I,J]

        return Mdense
    # ---------------------------------------


    # ---------------------------------------
    def genSpMatrix(n,m):
        from scipy import sparse

        randint = np.random.randint
        randv = np.random.random

        nnz = int(np.sqrt(n*m))

        I = np.array([randint(0,n) for i in range(0,nnz)])
        J = np.array([randint(0,m) for j in range(0,nnz)])
        V = randv(nnz)

        return sparse.coo_matrix((V,(I,J)),shape=(n,m)).tocsr()
    # ---------------------------------------

    # ---------------------------------------
    def genSpVector(n):
        from scipy import sparse
        import numpy as np

        randv = np.random.random

        m = 1
        nnz = n

        I = np.array(list(range(0,n)))
        J = np.array([0]*nnz)
        V = randv(nnz)

        return sparse.coo_matrix((V,(I,J)),shape=(n,m)).tocsr()
    # ---------------------------------------

    n = 100
    m = 100

    Mx  = genSpMatrix(n,m)
    My  = genSpMatrix(n,m)
    S   = genSpMatrix(n,m)
    Dx  = genSpMatrix(n,m)
    Dy  = genSpMatrix(n,m)
    Rx  = genSpMatrix(n,m)
    Ry  = genSpMatrix(n,m)

    X   = genSpVector(n)
    Y   = genSpVector(n)

    from scipy.sparse import eye
    from scipy.io import mmwrite

    Zxy_ = 0. * eye(Mx.shape[0], My.shape[1]) ; Zxy = Zxy_.tocsr() # zero matrix
    Zyx_ = 0. * eye(My.shape[0], Mx.shape[1]) ; Zyx = Zyx_.tocsr() # zero matrix
    I_   = 1. * eye(1, 1) ; I = I_.tocsr()


    # ...
    matrices1 = [ [-S, -0.5*Rx, 0.5*Ry] \
                , [Dx,     -Mx,    Zxy] \
                , [Dy,     Zyx,    -My] ]


    matrices2 = [ [Dx,     -Mx,    Zxy] \
                , [Dy,     Zyx,    -My] ]

    matrices3 = [  [S, Rx] \
                , [Dx, Mx] ]

    matrices4 = [ [Dx, Mx] ]

    matrices5 = [ [Dx] , [Dy] ]

    Yt = Y.transpose().tocsr()
    matrices6 = [ [S , X] \
                , [Yt, I] ]

    matrices7 = [ [-S, -0.5*Rx, 0.5*Ry] \
                , [Dx,     -Mx,    None] \
                , [Dy,     None,    -My] ]
    # ...

    list_matrices = []
    list_matrices += [matrices1]
    list_matrices += [matrices2]
    list_matrices += [matrices3]
    list_matrices += [matrices4]
    list_matrices += [matrices5]
    list_matrices += [matrices6]
    list_matrices += [matrices7]
    for i,matrices in enumerate(list_matrices):
        print("============= Test "+ str(i+1) +" ===========")

        A = BlockMatrix(matrices)
        A.assembly()

        Adense = assembly_dense(A.matrices)

        print("Assembling error = ", np.sum(A.todense()-Adense))

# -*- coding: UTF-8 -*-
"""
This module is intend to solve the matrix equation

sum_i=1^r  Ax_i X Ay_i = F

Where   Ay_i are circulant matrices
        Ax_i is in general are band matrices
"""

import numpy as np
from scipy.linalg import circulant, inv
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import gmres, splu
from scipy.sparse import kron
from scipy.io import mmwrite, mmread
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator


# -----------------------
pi =  np.pi
cos = np.cos
sin = np.sin
# -----------------------

# ...
def CMPLX(x,y):
    return x + y * 1j
# ...

# ...
def genTestMatrices(r, nx, ny, p, EXPORT=False, IMPORT=False):
    list_Ax = [] ; list_Ay = []

    # ... Define non singular diagonal matrices for the x-direction
    shift = 0
    for i in range(0,r):
        if IMPORT:
            Ax = mmread("figa/Ax"+str(i)+".mtx")
        else:
#           a = np.random.random(nx)
            a = np.ones(nx)
            Ax = diags(a,shift)
        Ax = csr_matrix(Ax)
        if EXPORT:
            mmwrite("figa/Ax"+str(i)+".mtx", Ax)

        list_Ax.append(Ax)
    # ...

    # ... Define circulant matrices for the x-direction
    for i in range(0,r):
        if IMPORT:
            Ay = mmread("figa/Ay"+str(i)+".mtx")
        else:
            ay = np.zeros(ny)
            ay[:2*p+1] = np.random.random(2*p+1)
            Ay = circulant(ay)
        Ay = csr_matrix(Ay)
        if EXPORT:
            mmwrite("figa/Ay"+str(i)+".mtx", Ay)

        list_Ay.append(Ay)
    # ...

    return list_Ax, list_Ay
# ...

# ...
def computeEigenValues(list_Ay, cmplx=True):
    # ...
    def computeEigenVal(A):
        dtype = np.double
        if cmplx:
            dtype = np.complex
        n,m = A.shape
        a = np.zeros(n)
        for i in range(0,n):
            a[i] = A[0,i]
        eigenA = np.zeros(n, dtype=dtype)
        for k in range(0,n):
            ck = 2*pi*k/n
            for j in range(0,n):
                if cmplx:
                    eigenA[k] += a[j] * CMPLX( cos(ck*j) , sin(ck*j) )
                else:
                    eigenA[k] += a[j] * cos(ck*j)
        return eigenA
    # ...
    list_eigenAy = []
    for Ay in list_Ay:
        eigenAy = computeEigenVal(Ay)
        list_eigenAy.append(eigenAy)
    return list_eigenAy
# ...

# ...
def AssembleColumnMatrix(j, nx, ny, list_Ax, list_eigenAy):
    """
    j must be in range(0,ny)
    """
    Sp = np.zeros((nx,nx))
    for Ax,eigenAy in zip(list_Ax,list_eigenAy):
        Sp = Sp + eigenAy[j] * Ax.todense()
    return csr_matrix(Sp)
# ...

# ...
def solveSp(Sp, b):
    x = gmres(Sp,b)[0]
    return x
# ...

# ...
def rsolve(list_Ax, list_eigenAy, F):
    fft = np.fft.rfft
    ifft = np.fft.irfft

    # ...
    nx,ny = F.shape
    n = nx ; m = ny
    mmax = m/2 -1
    x = F.transpose()
    _F = np.zeros((m, n))
    U = np.zeros_like(_F)
    # ...

    # ...
    y = np.zeros((m/2 + 1, n), dtype=np.complex)
    for j in range(0, n):
        x1d = x[:,j]
        y1d = fft(x1d)
        y[:,j] = y1d
    # ...

    # ... if ny is even
    for j in range(0,n):
        _F[0,j] = y[0,j].real
        for i in range(1,mmax+1):
            z = y[i,j]
            _F[2*i-1,j] = z.real
            _F[2*i  ,j] = z.imag
        _F[m-1,j] = y[m/2,j].real
    # ...

    # ...

    # ... treatment of the 0-mode
    f1d = _F[0, :]
    Sp = AssembleColumnMatrix(0, nx, ny, list_Ax, list_eigenAy)
    u1d = solveSp(Sp, f1d)
    U[0, :] = u1d

    for j in range(1, mmax+1):
        Sp = AssembleColumnMatrix(j, nx, ny, list_Ax, list_eigenAy)

        # ... treatment of the mode 2j-1
        f1d = _F[2*j-1, :]
        u1d = solveSp(Sp, f1d)
        U[2*j-1, :] = u1d

        # ... treatment of the mode 2j
        f1d = _F[2*j, :]
        u1d = solveSp(Sp, f1d)
        U[2*j, :] = u1d

    # ... treatment of the last mode
    f1d = _F[m-1, :]
    Sp = AssembleColumnMatrix(mmax+1, nx, ny, list_Ax, list_eigenAy)
    u1d = solveSp(Sp, f1d)
    U[m-1, :] = u1d
    # ...

    # ... if ny is even
    y = np.zeros_like(y)
    for j in range(0,n):
        y[0, j] = CMPLX(U[0, j], 0.0)

        for i in range(1, mmax+1):
            y[i, j] = CMPLX ( U[2*i - 1, j] , U[2*i, j] )

        y[m/2, j] = CMPLX ( U[m-1, j] , 0.0 )
    # ...

    # ...
    x = np.zeros_like(x)
    for j in range(0, n):
        y1d = y[:,j]
        x1d = ifft(y1d)
        x[:,j] = x1d
    # ...

    # ...
    X = x.transpose()
    #print X
    # ...

    return X
# ...

# ...
def csolve(list_Ax, list_eigenAy, F, EXPORT=False, list_opSj=None):
    fft = np.fft.fft
    ifft = np.fft.ifft

    X  = np.zeros_like(F)
    Yp = np.zeros_like(F, dtype=np.complex)
    Xp = np.zeros_like(F, dtype=np.complex)

    # ...
    for i in range(0, nx):
        # ... extract the i^th line as a vector
        y = F[i,:]
        # ... move to the commun basis using FFT
        yp = fft(y)
        Yp[i,:] = yp
    # ...

    # ...
    for j in range(0, ny):
        if list_opSj is None:
            # ... assemble the 1D matrix
            Sj = AssembleColumnMatrix(j, nx, ny, list_Ax, list_eigenAy)
        if EXPORT:
            mmwrite("figa/S"+str(j)+".mtx", Sj)
        # ... extract the j^th column as a vector
        yp = Yp[:,j]
        # ... solve the 1D linear system in the commun basis
        if list_opSj is None:
            xp = gmres(Sj,yp)[0]
        else:
            opSj = list_opSj[j]
            xp = opSj.solve(yp)
        Xp[:,j] = xp
    # ...

    # ...
    for i in range(0, nx):
        xp = Xp[i,:]
        # ... come back to the real space
        x = ifft(xp)
        # ... ... make sur that it is real
        x = x.real
        # ... update the global matrix
        X[i,:] = x
    # ...

    return X
# ...


# ...
def verification(list_Ax, list_Ay, X, F):
    _F = np.zeros_like(X)
    for Ax,Ay in zip(list_Ax, list_Ay):
        _F += Ax * X * Ay.transpose()
#    print "F   ", F
#    print "_F  ", _F
    print((np.allclose(F, _F)))
#    assert(np.allclose(F, _F))
# ...

# ...
def constructGlobalSystem(list_Ax, list_Ay):
    # ...
    list_eigenAy = computeEigenValues(list_Ay)
    # ...

    # ...
    Ax0 = list_Ax[0]
    Ay0  = list_Ay[0]
    S = kron(Ay0, Ax0)
    r = len(list_Ax)
    for i in range(1, r):
        Ax = list_Ax[i]
        Ay  = list_Ay[i]
        S  = S + kron(Ay, Ax)
    return S
    # ...
# ...

# ...
class nearestCirculant(object):
    """
    this class constructs a list of circulant matrices that approche a given
    list of matrices A by minimizing the Frobenius norm
    """
    def __init__(self, list_A, cost=0):
        self.list_A = list_A
        self.method = method

        norm = lambda M: np.linalg.norm(M, 'fro')

        # ...
        def cost0(M, c):
            C = circulant(c)
            nr = norm(M-C)
            return nr
        # ...

        # ...
        def cost1(M, c):
            n,m = M.shape
            C = circulant(c)
            invC = inv(C)
            I = np.eye(n)
            nr = norm(I-invC*M)
            return nr
        # ...

        # ...
        def cost2(M, c):
            diag = M.diagonal()
            shift = 0
            D = diags(diag,shift)
            Z = M-D
            C = circulant(c)
            nr = norm(Z-C)
            return nr
        # ...

        self.cost0 = cost0
        self.cost1 = cost1
        self.cost2 = cost2
        self.cost = getattr(self, 'cost%d' % cost)

    def construct(self, method='BFGS', tol = 1.e-7):
        list_C = []
        for A in self.list_A:
            # ...
            if method is None:
                n,m = A.shape
                MD = A.todense()
                c = np.zeros(n)
                for k in range(0,n):
                    c1 =0.; c2=0.
                    for i in range(0,n-k):
                        c1 += MD[i,k+i]
                    for i in range(n-k,n):
                        c2 += MD[i,k+i-n]
                    c[k] = ( c1 + c2 ) / n
            else:
                cost = lambda c: self.cost(A,c)

                n,m = A.shape
                x0 = np.zeros(n)
                x0[0] = 1.
                res = minimize(  cost, x0 \
                               , method=method \
                               , options={'gtol': tol, 'disp': verbose})
                c = res.x
            # ...

            C = circulant(c)
            C = csr_matrix(C)
            list_C.append(C)

        return list_C
# ...

# ...
class circulantPrecond(object):
    def __init__(self, list_Ax, list_Ay \
                 , cost=0, method='BFGS' \
                 , tol = 1.e-7, verbose=False):

        # ... construct the nearest circulant matrices for list_Ay
        nearCirc = nearestCirculant(list_Ay, cost=cost)
        list_C  = nearCirc.construct(method=method, tol=tol)
        # ...

        self.list_C = list_C

        # ...
        self.list_eigenC        = computeEigenValues(list_C)
        # ...

        # ...
        n,m = list_Ax[0].shape ; nx = n
        n,m = list_Ay[0].shape ; ny = n
        self.n = [nx,ny]
        # ...

        # ...
        r = len(list_Ax)
        Ax0 = list_Ax[0]
        C0  = list_C[0]
        P = kron(C0, Ax0)
        for i in range(1, r):
            Ax = list_Ax[i]
            C  = list_C[i]
            P  = P + kron(C, Ax)
        self.P = P
        # ...

        # ...
        list_opSj = []
        for j in range(0, ny):
            # ... assemble the 1D matrix
            Sj = AssembleColumnMatrix(j, nx, ny, list_Ax, self.list_eigenC)
            opSj = splu(Sj.tocsc())
            list_opSj.append(opSj)
        self.list_opSj = list_opSj
        # ...

    def aspreconditioner(self):
        """Create a preconditioner

        Returns
        -------
        precond : LinearOperator
            Preconditioner suitable for the iterative solvers in defined in
            the scipy.sparse.linalg module (e.g. cg, gmres) and any other
            solver that uses the LinearOperator interface.  Refer to the
            LinearOperator documentation in scipy.sparse.linalg

        See Also
        --------
        scipy.sparse.linalg.LinearOperator

        Examples
        --------
        >>>

        """
        shape = self.P.shape
        dtype = self.P.dtype

        nx, ny = self.n

        self.i = 0
        def matvec(b):
            F = b.reshape((ny,nx))
            F = F.transpose()
            X = csolve(self.list_C, self.list_eigenC, F, list_opSj=self.list_opSj)
            x = X.transpose().reshape(nx*ny)
#            print ">> iteration ", self.i
            self.i += 1
            return x

        return LinearOperator(shape, matvec, dtype=dtype)
# ...

# ...
def testcase(r, nx, ny, p, EXPORT=False, IMPORT=False):
    # ...
    if IMPORT:
        F = np.genfromtxt("figa/F.txt")
        try:
            nx,ny = F.shape
        except:
            nx = 1
            ny, = F.shape
            _F = F
            F = np.zeros((nx,ny))
            F[0,:] = _F
    else:
        F = np.random.random((nx,ny))
        np.savetxt("figa/F.txt", F)
    # ...

    # ...
    list_Ax, list_Ay    = genTestMatrices(r, nx, ny, p \
                                          , EXPORT=EXPORT \
                                          , IMPORT=IMPORT)
    # ...

    return list_Ax, list_Ay, F
# ...

# ...
def testcase_poisson(scale=False):
    Mx = mmread("figa/Mx.mtx") ; Mx = Mx.tocsr()
    Sx = mmread("figa/Sx.mtx") ; Sx = Sx.tocsr()
    Kx = mmread("figa/Kx.mtx") ; Kx = Kx.tocsr()
    KTx = Kx.transpose().tocsr()

    My = mmread("figa/My.mtx") ; My = My.tocsr()
    Sy = mmread("figa/Sy.mtx") ; Sy = Sy.tocsr()
    Ky = mmread("figa/Ky.mtx") ; Ky = Ky.tocsr()
    KTy = Ky.transpose().tocsr()

#    # ...
#    list_Ax = [Mx, Sx,  Kx, KTx]
#    list_A  = [Sy, My, KTy,  Ky]
#    # ...

#    # ...
#    Kmx =   np.sqrt(2) * (Kx+KTx)
#    Kjx = - np.sqrt(2) * (Kx-KTx)
#
#    Kmy =   np.sqrt(2) * (Ky+KTy)
#    Kjy =   np.sqrt(2) * (Ky-KTy)
#
#    list_Ax = [Mx, Sx, Kmx, Kjx]
#    list_A  = [Sy, My, Kmy, Kjy]
#    # ...

#    # ...
#    list_Ax = [ Kx, KTx, Sx]
#    list_A  = [KTy,  Ky, My]
#    # ...

    # ...
    list_Ax = [Mx, Sx]
    list_A  = [Sy, My]
    # ...

    if scale:
        print("MUST IMPROVED: WE HAVE TO MULTIPLY BY ONE MATRIX FOR ALL MATRICES")
        shift = 0
        list_Ay = []
        for A in list_A:
            diag = 1./A.diagonal()
            D = diags(diag, shift).tocsr()
            Ay = A * D
            Ay.tocsr()
            list_Ay.append(Ay)
    else:
        list_Ay = list_A

    n,m = Mx.shape ; nx = n
    n,m = My.shape ; ny = n
#    F = np.random.random((nx,ny))
    F = np.ones((nx,ny))

    return list_Ax, list_Ay, F
# ...

# ---------------------------------------------------------------
if __name__=="__main__":
    from time import time
    # -------------------------
#    nx = 512 ; ny = 512
#    nx = 256 ; ny = 256
#    nx = 128 ; ny = 128
#    nx =  64 ; ny =  64
    nx =  32 ; ny =  32
#    nx =  16 ; ny =  16

    r       = 4
    p       = 3

#    EXPORT = True
    EXPORT  = False

    IMPORT  = False
#    IMPORT = True

    method  = None
    cost    = 0
#    method  = 'BFGS'
    tol     = 1.e-7
#    verbose = True
    verbose = False

#    scale   = True
    scale   = False

#    CIRCULANT = True
    CIRCULANT = False
    # -------------------------

    # ...
    if CIRCULANT:
        list_Ax, list_Ay, F = testcase(r, nx, ny, p, EXPORT=False, IMPORT=False)
    else:
        list_Ax, list_Ay, F = testcase_poisson(scale=scale)
#        n,m = list_Ax[0].shape
#        r = len(list_Ax)
#        list_Ax = []
#        for i in range(0,r):
##            diag = np.random.random(n)
#            diag = np.ones(n)
#            shift = 0
#            A = diags(diag, shift)
#            list_Ax.append(A)

#        _list_Ax = list_Ax[:3]
#        _list_Ay = list_Ay[:3]

        _list_Ax = list_Ax[:2]
        _list_Ay = list_Ay[:2]
        PrecConstruct = circulantPrecond(_list_Ax, _list_Ay \
                                         , cost=cost, method=method \
                                         , tol=tol, verbose=verbose)
        mmwrite('figa/P.mtx', PrecConstruct.P)

#        mmwrite('figa/C_Sy.mtx', PrecConstruct.list_C[0])
#        mmwrite('figa/C_My.mtx', PrecConstruct.list_C[1])
#        mmwrite('figa/C_Kmy.mtx', PrecConstruct.list_C[2])
#        mmwrite('figa/Kmy.mtx', list_Ay[2])

#        mmwrite('figa/C_KTy.mtx', PrecConstruct.list_C[2])
#        mmwrite('figa/C_Ky.mtx' , PrecConstruct.list_C[3])

#        mmwrite('figa/C_KTy.mtx', PrecConstruct.list_C[0])
#        mmwrite('figa/C_Ky.mtx' , PrecConstruct.list_C[1])
#        mmwrite('figa/C_My.mtx' , PrecConstruct.list_C[2])

        mmwrite('figa/C_Sy.mtx', PrecConstruct.list_C[0])
        mmwrite('figa/C_My.mtx', PrecConstruct.list_C[1])

        Precond = PrecConstruct.aspreconditioner()
    # ...

    # ...
    n,m = list_Ax[0].shape ; nx = n
    n,m = list_Ay[0].shape ; ny = n
    # ...

    # ...
    S = constructGlobalSystem(list_Ax, list_Ay)
    mmwrite('figa/S.mtx', S)
    # ...

    # ...
    print("=============================")
    print("  nx, ny  ", nx, ny)
    print("  size    ", S.shape)
    print("  nnz     ", S.nnz)
    print("=============================")
    # ...

#    import sys ; sys.exit(0)

    # ...
    print("=============================")
    print(">>> using  the global system")
    y = F.transpose().reshape(nx*ny)
    tb = time()
    Xg,it = gmres(S, y)
    Xg = Xg.reshape((ny,nx))
    Xg = Xg.transpose()
    te = time()
    print("Elapsed time ", te-tb)
    # ...

    # ...
    if CIRCULANT:
        print("=============================")
        print(">>> using circulant fast solver")
        list_eigenAy  = computeEigenValues(list_Ay)
        tb = time()
        X = csolve(list_Ax, list_eigenAy, F)
        te = time()
        print("Elapsed time ", te-tb)
        print("Internal verification ")
        verification(list_Ax, list_Ay, X, F)
    else:
        print("=============================")
        print(">>> using circulant preconditioner solver")
        tb = time()
        y = F.transpose().reshape(nx*ny)
        x,it = gmres(S, y, M=Precond)
        X = x.reshape((ny,nx))
        X = X.transpose()
        te = time()
        print("Elapsed time ", te-tb)
    # ...

    # ...
    print("=============================")
    print("Is everything OK?")
    print(np.allclose(Xg,X, rtol=1e-07) \
            , " with error ", np.linalg.norm(Xg-X)/np.linalg.norm(X))


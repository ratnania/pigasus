# -*- coding: UTF-8 -*-
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from caid.core.bspline import bsp
from sys import exit as STOP

# ----------------------------
xb = 0. ; xe = 1.

#p = 2
p = 3

#N = 63
N = 31

#N = 128
#N = 64
#N = 32
#N = 16
#N = 8

Nnew = N+p+1

PLOT = True
PLOT = False
NPTS = 200
# ----------------------------

# ----------------------------
T       = list(range(-p,0)) + list(range(0,N)) + list(range(N,N+p+1))
T       = np.asarray(T, dtype=np.double) / N
Interv  = list(range(0,N+1))
Interv  = np.asarray(Interv, dtype=np.double) / N
# ----------------------------

# ----------------------------
def genGrid(Interv, p, type="legendre"):
    from pigasus.fem.quadratures import quadratures
    qd = quadratures()
    list_xgl, list_wgl = qd.generate(Interv, p, as_type=type)
    xgl = [] ; wgl = []
    for x,w in zip(list_xgl, list_wgl):
        xgl += list(x)
        wgl += list(w)
    xgl = np.asarray(xgl)
    wgl = np.asarray(wgl)

    return xgl, wgl
# ----------------------------

# ----------------------------
def plotBasis():
    X = np.linspace(xb,xe,NPTS)
    M = assembleM(X, T, p, N)
    list_i = list(range(0,N))
    for i in list_i:
        B = M[:,i]
        plt.plot(X, B)
    plt.show()
# ----------------------------

# ------------------------------
def findSpan(x, tx, px):
    leftmkx = bsp.FindSpan(px,tx,x) - px

    return leftmkx
# ------------------------------

# ------------------------------
def evalSplines(x, tx, px, deriv=0):
#    Nx = np.zeros(px+1)
#    Nx[:] = bsp.EvalBasisFuns(px,tx,x)
    Y = bsp.EvalBasisFunsDers(px,tx,x,deriv)
    return Y[deriv,:]
# ------------------------------

# ------------------------------
def assembleM(X, tx, px, nx, deriv=0):
    """
    """
    n, = X.shape

    # initialize the matrix M of size (nx , n )
    M = np.zeros((n, nx))

    for i in range(0,n):
#        print ">>>> i ", i
        I = i
        x = X[i]
        Nx      = evalSplines(x, tx, px, deriv=deriv)
        leftmkx = findSpan   (x, tx, px)
#        print '>> x  ', x
#        print 'Basis ', Nx
        for ip in range(0, px+1):
            J = ip + leftmkx
            if J >= nx:
                J -= nx
#            print "ip ", ip
#            print "I,J ", I,J

            M[I,J] += Nx[ip]

    return M
# ------------------------------

# ------------------------------
def assembleR(X, W, fct):
    """
    """
    M = W * fct(X)
    return M
# ------------------------------

# ------------------------------
def assemble(A, R ,B):
    """
    """
    n,nx = A.shape

    M = np.zeros((nx,nx))
    for j in range(0, nx):
        for i in range(0, nx):
            M[i,j] += np.sum( A[:,i] * R[:] * B[:,j] )

    return M
# ------------------------------

if __name__ == "__main__":
    from pigasus.gallery.basicPDE import *
    from caid.cad_geometry import cad_geometry
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix

    #-----------------------------------
    t = 0.5
#    nx = 2 ; px = 1
    nx = N-(p-1) ; px = p
    nlevel = 1
    #-----------------------------------

    geo = cad_geometry("iter_inner.xml")
    nrb = geo[0]

    #-----------------------------------
#    # ...
#    def a_ext(s):
#        S = nrb(s)
#        return S[:,0]
#    # ...
#    # ...
#    def b_ext(s):
#        S = nrb(s)
#        return S[:,1]
#    # ...
#    # ...
#    def ads_ext(s):
#        S = nrb.evaluate_deriv(s, nderiv=1)
#        return S[1,:,0]
#    # ...
#    # ...
#    def bds_ext(s):
#        S = nrb.evaluate_deriv(s, nderiv=1)
#        return S[1,:,1]
#    # ...
    # ...
    def a_ext(s):
        return np.cos(2*np.pi*s)
    # ...
    # ...
    def b_ext(s):
        return np.sin(2*np.pi*s)
    # ...
    # ...
    def ads_ext(s):
        return -2*np.pi*np.sin(2*np.pi*s)
    # ...
    # ...
    def bds_ext(s):
        return -2*np.pi*np.cos(2*np.pi*s)
    # ...
    # ...
    def Lambda11(s):
        ae = a_ext(s)
        be = b_ext(s)
        ae_ds = ads_ext(s)
        be_ds = bds_ext(s)

        G  = ae_ds**2   + be_ds**2
        G /= ae * be_ds - be * ae_ds

        return G
    # ...
    # ...
    def Lambda12(s):
        ae = a_ext(s)
        be = b_ext(s)
        ae_ds = ads_ext(s)
        be_ds = bds_ext(s)

        G  = ae * ae_ds + be * be_ds
        G /= ae * be_ds - be * ae_ds

        return G
    # ...
    # ...
    def Lambda22(s):
        ae = a_ext(s)
        be = b_ext(s)
        ae_ds = ads_ext(s)
        be_ds = bds_ext(s)

        G  = ae**2      + be**2
        G /= ae * be_ds - be * ae_ds

        return G
    # ...
#    s = np.linspace(0.,1.,100)
#    ae = a_ext(s)
#    be = b_ext(s)
#    ae_ds = ads_ext(s)
#    be_ds = bds_ext(s)
    # ...
    if PLOT:
        plotBasis()
        plt.plot(ae,be)
        plt.show()
    # ...
    #-----------------------------------

    #-----------------------------------
    # ... Computing Quadrature Points and Weights
    xgl, wgl = genGrid(Interv, p, type="legendre")
    # ...
    #-----------------------------------

    #-----------------------------------
    # ...
    A  = assembleM(xgl, T, p, N, deriv=0)
    R  = assembleR(xgl, wgl, Lambda11)
    B  = assembleM(xgl, T, p, N, deriv=0)

    M  = assemble(A, R ,B)
    # ...

    # ...
    A  = assembleM(xgl, T, p, N, deriv=1)
    R  = assembleR(xgl, wgl, Lambda22)
    B  = assembleM(xgl, T, p, N, deriv=1)

    S  = assemble(A, R ,B)
    # ...

    # ...
    A  = assembleM(xgl, T, p, N, deriv=1)
    R  = assembleR(xgl, wgl, Lambda12)
    B  = assembleM(xgl, T, p, N, deriv=0)

    K  = assemble(A, R ,B)
#    KT = assemble(B, R ,A)

    My = M ; Ky = K ; KTy = K.transpose() ; Sy = S
    np.savetxt("figa/My.txt" , My)
    np.savetxt("figa/Sy.txt" , Sy)
    np.savetxt("figa/Ky.txt" , Ky)
    mmwrite("figa/My.mtx" , csr_matrix(My))
    mmwrite("figa/Sy.mtx" , csr_matrix(Sy))
    mmwrite("figa/Ky.mtx" , csr_matrix(Ky))
    # ...
    #-----------------------------------

    #-----------------------------------
    # ...
    def Gamma11(s):
        G  = t**2 * (1-s)**2 + 2*t*s*(1-s) + s
        G /= t*(1-t)*(1-s) + (1-t)*s

        return [G]
    # ...
    # ...
    def Gamma12(s):
        G  = -1.

        return [G]
    # ...
    # ...
    def Gamma22(s):
        G  = (1-t)**2
        G /= t*(1-t)*(1-s) + (1-t)*s

        return [G]
    # ...
    #-----------------------------------

    #-----------------------------------
    # ...
    tc = {}
    tc['A']             = Gamma11
    tc['u']             = lambda x : [0.]
    tc['f']             = lambda x : [0.]
    tc['Dirichlet']     = [[0,1]]
    tc11 = tc
    # ...
    # ...
    tc = {}
    tc['v']             = Gamma12
    tc['u']             = lambda x : [0.]
    tc['f']             = lambda x : [0.]
    tc['Dirichlet']     = [[0,1]]
    tc12 = tc
    # ...
    # ...
    tc = {}
    tc['b']             = Gamma22
    tc['u']             = lambda x : [0.]
    tc['f']             = lambda x : [0.]
    tc['Dirichlet']     = [[0,1]]
    tc22 = tc
    # ...
    #-----------------------------------

    #-----------------------------------
    # ...
    from caid.cad_geometry import line as domain
    geo = domain(n=[nx],p=[px])
    # ...
    PDE11 = basicPDE(geometry=geo, testcase=tc11)
    PDE12 = basicPDE(geometry=geo, testcase=tc12)
    PDE22 = basicPDE(geometry=geo, testcase=tc22)

    PDE11.assembly()
    PDE12.assembly()
    PDE22.assembly()

    M = PDE22.system.get().todense()
    K = PDE12.system.get().todense()
    S = PDE11.system.get().todense()

    Mx = M ; Kx = K ; KTx = K.transpose() ; Sx = S
    np.savetxt("figa/Mx.txt" , Mx)
    np.savetxt("figa/Sx.txt" , Sx)
    np.savetxt("figa/Kx.txt" , Kx)
    mmwrite("figa/Mx.mtx" , csr_matrix(Mx))
    mmwrite("figa/Sx.mtx" , csr_matrix(Sx))
    mmwrite("figa/Kx.mtx" , csr_matrix(Kx))
    #-----------------------------------

    import sys ; sys.exit(0)

    #-----------------------------------
    import pywt
    import numpy as np

    Sx = np.genfromtxt("figa/Sx.txt")
    Mx = np.genfromtxt("figa/Mx.txt")
    Kx = np.genfromtxt("figa/Kx.txt")

    Sy = np.genfromtxt("figa/Sy.txt")
    My = np.genfromtxt("figa/My.txt")
    Ky = np.genfromtxt("figa/Ky.txt")

    KTx = Kx.transpose()
    KTy = Ky.transpose()

    list_Mx = [Sx, Kx , KTx, Mx]
    list_My = [My, KTy, Ky , Sy]


    # ...
    list_coeffsx    = [] ; list_coeffsy = []
    list_Ux         = [] ; list_Uy      = []
    nM = len(list_Mx)
    for (A,B) in zip(list_Mx, list_My):
        data = A
#        coeffs = pywt.dwt2(data, 'haar')
#        rdata = pywt.idwt2(coeffs, 'haar')
#        assert(np.allclose(rdata,data))

        coeffs = pywt.wavedec2(data, 'db1', level=nlevel)
        rdata = pywt.waverec2(coeffs, 'db1')
        assert(np.allclose(rdata,data))

        list_coeffsx.append(coeffs)
        list_Ux.append(coeffs[0])

        data = B
#        coeffs = pywt.dwt2(data, 'haar')
#        rdata = pywt.idwt2(coeffs, 'haar')
#        assert(np.allclose(rdata,data))

        coeffs = pywt.wavedec2(data, 'db1', level=nlevel)
        rdata = pywt.waverec2(coeffs, 'db1')
        assert(np.allclose(rdata,data))

        list_coeffsy.append(coeffs)
        list_Uy.append(coeffs[0])
    # ...

    # ...
    #-----------------------------------

    #-----------------------------------
    from scipy.linalg import kron
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import cg

    i = 0
    Ax = list_Mx[i]
    Ay = list_My[i]
    M = kron(Ax,Ay)
    for i in range(1, len(list_Mx)):
        Ax = list_Mx[i]
        Ay = list_My[i]
        M += kron(Ax,Ay)

    data = M
    coeffs = pywt.wavedec2(data, 'db1', level=2*nlevel)
    rdata = pywt.waverec2(coeffs, 'db1')
    assert(np.allclose(rdata,data))
    cM = coeffs[0]

    M = csr_matrix(M)

    n,m = M.shape
    b = np.ones(n)

    x = cg(M, b)[0]
    data    = x
    coeffs  = pywt.wavedec(data, 'db1', level=2*nlevel)
    rdata   = pywt.waverec(coeffs, 'db1')
    assert(np.allclose(rdata,data))
#    print "x            ", x
#    print "compressed x ", coeffs[0]
    #-----------------------------------

    #-----------------------------------
    i = 0
    Ux = list_Ux[i]
    Uy = list_Uy[i]
    U = kron(Ux,Uy)
    for i in range(1, len(list_Mx)):
        Ux = list_Ux[i]
        Uy = list_Uy[i]
        U += kron(Ux,Uy)

#    print "Norm ", np.linalg.norm(cM, 'fro')
#    print "Norm ", np.linalg.norm(U, 'fro')
#    print np.allclose(cM, U, rtol=1e-03)
#    print "---"
#    print cM[0,:]
#    print "---"
#    print U[0,:]
#    print "---"

    Mat = csr_matrix(U)

    data    = b
    coeffs  = pywt.wavedec(data, 'db1', level=2*nlevel)
    rdata   = pywt.waverec(coeffs, 'db1')
    assert(np.allclose(rdata,data))

    u = coeffs[0]

#    print Mat.shape
#    print u.shape
    y = cg(Mat, u)[0]
#    print "compressed solution ", y
    #-----------------------------------

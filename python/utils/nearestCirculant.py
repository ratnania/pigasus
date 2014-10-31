# -*- coding: UTF-8 -*-
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import kron, circulant, inv
from scipy.sparse.linalg import cg
from scipy.sparse import csr_matrix, diags
from scipy.io import mmwrite, mmread

# ----------------
#MINIMIZE = True
MINIMIZE = False

COST = 0

#method = 'CG'
method = 'BFGS'

tol = 1.e-7
verbose = True

MAXITER = 100
FREQ    = 10
list_iters = list(range(0,MAXITER,FREQ))[1:]
# ----------------

norm = lambda M: np.linalg.norm(M, 'fro')

# ...
M = mmread("figa/My.mtx")
M = M.tocsr()
# ...

# ...
#from scipy.linalg import toeplitz
#import numpy as np
#p = 3
#n = 128
#a=np.zeros(n) ; b=np.zeros(n)
#a[:p+1] = np.random.random(p+1)
#a[n-p:] = np.random.random(p)
#b[:p+1] = np.random.random(p+1)
#b[n-p:] = np.random.random(p)
#a[0] = 1.
#b[0] = a[0]
#M = toeplitz(a,b)
# ...

# ...
if MINIMIZE:
    M = csr_matrix(M)
    diag = M.diagonal()
    shift = 0
    D = diags(diag,shift)
    Z = M-D
    n,m = M.shape
    # ...

    # ...
    def cost0(c):
        C = circulant(c)
        nr = norm(M-C)
        return nr
    # ...

    # ...
    def cost1(c):
        C = circulant(c)
        invC = inv(C)
        I = np.eye(n)
        nr = norm(I-invC*M)
        return nr
    # ...

    # ...
    def cost2(c):
        C = circulant(c)
        nr = norm(Z-C)
        return nr
    # ...

    # ...
    if COST == 0:
        cost = cost0
    if COST == 1:
        cost = cost1
    if COST == 2:
        cost = cost2
    # ...

    # ...
    x0 = np.zeros(n)
    x0[0] = 1.
    res = minimize(  cost, x0 \
                   , method=method \
                   , options={'gtol': tol, 'disp': verbose})
    c = res.x
    # ...
else:
    # ...
    n,m = M.shape
    MD = M.todense()
    c = np.zeros(n)
    for k in range(0,n):
        c1 =0.; c2=0.
        for i in range(0,n-k):
            c1 += MD[i,k+i]
        for i in range(n-k,n):
            c2 += MD[i,k+i-n]
        c[k] = ( c1 + c2 ) / n
    # ...

# ...
C = circulant(c)
invC = inv(C)
C = csr_matrix(C)
# ...

# ...
if COST in [0,1]:
    P = invC
if COST in [2]:
    _P = D + C
    _P = _P.todense()
    P = inv(_P)

#    idiag = np.array([1./d for d in diag])
#    shift = 0
#    invD = diags(idiag,shift)
#    P = invD + invC
# ...


# ...
b = np.ones(n)

print("----- stand-alone     cg ----")
for maxiter in list_iters:
    x,niter = cg(M, b, maxiter=maxiter)
    print("Error after ", maxiter, " iterations : ", np.linalg.norm(M*x-b))

print("----- preconditionned cg ----")
for maxiter in list_iters:
    x,niter = cg(M, b, M=P, maxiter=maxiter)
    print("Error after ", maxiter, " iterations : ", np.linalg.norm(M*x-b))
# ...


mmwrite("/home/ratnani/M.mtx", M)
mmwrite("/home/ratnani/P.mtx", C)

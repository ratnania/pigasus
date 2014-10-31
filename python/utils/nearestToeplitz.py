# -*- coding: UTF-8 -*-
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import kron, toeplitz, inv
from scipy.sparse.linalg import cg
from scipy.sparse import csr_matrix, diags

# ----------------
COST = 2

#method = 'CG'
method = 'BFGS'

tol = 1.e-7
verbose = True

MAXITER = 150
FREQ    = 10
list_iters = list(range(0,MAXITER,FREQ))[1:]
# ----------------

norm = lambda M: np.linalg.norm(M, 'fro')
#norm = lambda M: np.linalg.norm(M, 1)
#norm = lambda M: np.linalg.norm(M, 2)
#norm = lambda M: np.linalg.norm(M, -1)

M = np.genfromtxt("My.txt")
M = csr_matrix(M)
diag = M.diagonal()
shift = 0
D = diags(diag,shift)
Z = M-D
n,m = M.shape
# ...

# ...
def cost0(c):
    a = c[:n] ; b = c[n:]
    T = toeplitz(a,b)
    nr = norm(M-T)
    return nr
# ...

# ...
def cost1(c):
    a = c[:n] ; b = c[n:]
    T = toeplitz(a,b)
    invT = inv(T)
    I = np.eye(n)
    nr = norm(I-invT*M)
    return nr
# ...

# ...
def cost2(c):
    a = c[:n] ; b = c[n:]
    T = toeplitz(a,b)
    nr = norm(Z-T)
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


x0 = np.zeros(2*n)
x0[0] = 1.
res = minimize(  cost, x0 \
               , method=method \
               , options={'gtol': tol, 'disp': verbose})
c = res.x
a = c[:n] ; b = c[n:]
T = toeplitz(a,b)
invT = inv(T)
T = csr_matrix(T)
# ...

# ...
if COST in [0,1]:
    P = invT
if COST in [2]:
    _P = D + T
    _P = _P.todense()
    P = inv(_P)

#    idiag = np.array([1./d for d in diag])
#    shift = 0
#    invD = diags(idiag,shift)
#    P = invD + invT
# ...

# ...
b = np.ones(n)

print("###############################################")
print("        Solving ", n, "x" , n ," Problem")
print("###############################################")

print("----- stand-alone     cg ----")
for maxiter in list_iters:
    x,niter = cg(M, b, maxiter=maxiter)
    print("Error after ", maxiter, " iterations : ", np.linalg.norm(M*x-b))

print("----- preconditionned cg ----")
for maxiter in list_iters:
    x,niter = cg(M, b, M=P, maxiter=maxiter)
    print("Error after ", maxiter, " iterations : ", np.linalg.norm(M*x-b))
# ...


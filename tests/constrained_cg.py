from scipy.io import mmread
from scipy.sparse.linalg import cg, LinearOperator
from scipy.linalg import det
import numpy as np
try:
    from pyamg.util import condest
    withpyAMG = True
except:
    withpyAMG = False

A = mmread("A.mtx")
mean = np.genfromtxt("mean.txt")

if withpyAMG:
    print("Conditionning number ", condest(A))
print("Determinant          ", det(A.todense()))


n = A.shape[0]
def mv(v):
    r = A*v + np.dot(mean, v)
    return r

N = n
lin = LinearOperator((N,N), matvec=mv, dtype=np.float)

x = np.random.random(n)
b = A*x

x = np.zeros(N)
x1 = cg(A, b, tol=1.e-12)[0]
x2 = cg(lin, b, tol=1.e-12)[0]

print("|A*x1 - b|    ", np.linalg.norm(A*x1 - b))
print("mean(x1)      ", np.dot(mean, x1))
print("|L x2 - b|    ", np.linalg.norm(lin.matvec(x2) - b))
print("mean(x2)      ", np.dot(mean, x2))

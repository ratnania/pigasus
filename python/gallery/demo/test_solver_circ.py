import numpy as np
from figa_circ import *
from scipy.linalg import circulant, det
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu
from scipy.io import mmwrite, mmread

# -----------------------
pi =  np.pi
cos = np.cos

fft = np.fft.fft
ifft = np.fft.ifft
# -----------------------

# -------------------------
nx = 1 ; ny = 4
r = 1
p = 1
IMPORT = True
# -------------------------

# ...
if  IMPORT:
    A = mmread("figa/Ay0.mtx").tocsr()
    list_Ay = [A]
    list_Ax = []
    F = np.genfromtxt("figa/F.txt")
    b = F[:]
else:
    list_Ax, list_Ay    = genTestMatrices(r, nx, ny, p)
    b = np.random.random(ny)

list_eigenAy        = computeEigenValues(list_Ay, cmplx=True)
# ...

A = list_Ay[0]
eigA = list_eigenAy[0]

#A = A.todense()

# ...
print "======================="
print "Eigen values of A ", eigA
print "determinant  of A ", det(A.todense())
print "F  ", b
print "====================="
# ...

# ...
print "======================="
print ">>> Solving using FFT "
bp = fft(b)
xp = bp / eigA
x = ifft(xp)
x = x.real
print "x      ", x
print "Ax-b   ", np.linalg.norm(A.dot(x) - b)
# ...

# ...
print "======================="
print ">>> Solving using gmres "
x,it = gmres(A,b)
print "x      ", x
print "Ax-b   ", np.linalg.norm(A.dot(x) - b)
# ...

# ...
print "======================="
print ">>> Solving using splu "
opA = splu(A.tocsc())
x = opA.solve(b)
print "x      ", x
print "Ax-b   ", np.linalg.norm(A.dot(x) - b)
# ...

# ...
#print "======================="
#print ">>> Solving using invA "
#from scipy.linalg import inv
#invA = inv(A)
##from scipy.linalg import norm
##print A*invA
#x = invA * b
#print "x      ", x
#print "Ax-b   ", A.dot(x) - b
# ...



from scipy.sparse.linalg import cg, cgs, bicg, bicgstab, gmres, splu, spsolve
from pigasus.gallery.poisson import poisson
from igakit.cad_geometry import square as domain
import numpy as np
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

tol = 1.e-12
maxiter = 6000

#-----------------------------------
AllDirichlet = True

try:
    nx = int(sys.argv[1])
except:
    nx = 31

try:
    ny = int(sys.argv[2])
except:
    ny = 31

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2
#-----------------------------------

#-----------------------------------
geo = domain(n=[nx,ny],p=[px,py])
#-----------------------------------

# ...
try:
    bc_dirichlet
except NameError:
    bc_dirichlet = None
else:
    pass

try:
    bc_neumann
except NameError:
    bc_neumann = None
else:
    pass

try:
    AllDirichlet
except NameError:
    AllDirichlet = None
else:
    pass

try:
    Metric
except NameError:
    Metric = None
else:
    pass
# ...

PDE = poisson(geometry=geo, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann,
              AllDirichlet=AllDirichlet, metric=Metric)
PDE.assembly()
PDE.solve()

# getting scipy matrix
A = PDE.system.get()

b = np.ones(PDE.size)

print "Using cg."
x = cg(A, b, tol=tol, maxiter=maxiter)

print "Using cgs."
x = cgs(A, b, tol=tol, maxiter=maxiter)

print "Using bicg."
x = bicg(A, b, tol=tol, maxiter=maxiter)

print "Using bicgstab."
x = bicgstab(A, b, tol=tol, maxiter=maxiter)

print "Using gmres."
x = gmres(A, b, tol=tol, maxiter=maxiter)

print "Using splu."
op = splu(A.tocsc())
x = op.solve(b)

PDE.free()

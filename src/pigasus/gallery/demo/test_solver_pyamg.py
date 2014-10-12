from pigasus.gallery.poisson import poisson
from igakit.cad_geometry import square as domain
import numpy as np
import time
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

#-----------------------------------
tol_pyamg = 1.e-10
maxiter_pyamg = 200
accel       = 'gmres'
#accel       = 'cg'
maxiter_prec= 2
#-----------------------------------

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
A_scipy = PDE.system.get()

b = np.ones(PDE.size)

# ----------------------------------------------
import scipy
from pyamg.aggregation import smoothed_aggregation_solver

B = None                                # no near-null spaces guesses for SA

# Construct solver using AMG based on Smoothed Aggregation (SA) and display info
mls = smoothed_aggregation_solver(A_scipy, B=B)


# Solve Ax=b with no acceleration ('standalone' solver)
print "Using pyamg-standalone"
standalone_residuals = []
t_start = time.time()
x = mls.solve(b, tol=tol_pyamg, accel=None, maxiter=maxiter_pyamg, residuals=standalone_residuals)
t_end = time.time()
mls_elapsed = t_end - t_start
mls_err = standalone_residuals[-1]
mls_niter = len(standalone_residuals)
print "done."
standalone_residuals  = np.array(standalone_residuals)/standalone_residuals[0]
factor1 = standalone_residuals[-1]**(1.0/len(standalone_residuals))
standalone_final_err = np.linalg.norm(b-A_scipy.dot(x))


# Solve Ax=b with Conjugate Gradient (AMG as a preconditioner to CG)
print "Using pyamg-accelerated"
accelerated_residuals = []
t_start = time.time()
x = mls.solve(b, tol=tol_pyamg, accel=accel, maxiter=maxiter_pyamg, residuals=accelerated_residuals)
t_end = time.time()
mlsaccel_elapsed = t_end - t_start
mlsaccel_err = accelerated_residuals[-1]
mlsaccel_niter = len(accelerated_residuals)
print "done."
# Compute relative residuals
accelerated_residuals = np.array(accelerated_residuals)/accelerated_residuals[0]
# Compute (geometric) convergence factors
factor2 = accelerated_residuals[-1]**(1.0/len(accelerated_residuals))
accelerated_final_err = np.linalg.norm(b-A_scipy.dot(x))
# ----------------------------------------------

print "--------------- elapsed times ---------------"
print "mls          : ", mls_elapsed
print "mls-accel    : ", mlsaccel_elapsed
print "---------------------------------------------"

print "----------- PyAMG Multigrid solver ----------"
print mls
print "---------------------------------------------"

print "---------------------------------------------"
print "mls converges with error : ", mls_err, " after ", mls_niter, " iterationr"
print "final error : ", standalone_final_err
print "mlsaccel converges with error : ", mlsaccel_err, " after ", mlsaccel_niter, " iterations"
print "final error : ", accelerated_final_err
print "                     MG convergence factor: %g"%(factor1)
print "MG with CG acceleration convergence factor: %g"%(factor2)
print "---------------------------------------------"

PDE.free()

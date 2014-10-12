#! /usr/bin/python

import sys
import time
import numpy as np
from pigasus.gallery.monge_ampere import picard
from igakit.cad_geometry import square as domain
from matplotlib import pylab as plt
from pigasus.multigrid.agregation import *

sin = np.sin ; cos = np.cos ; exp = np.exp ; log = np.log ; sqrt = np.sqrt ; pi= np.pi ; sqrt = np.sqrt
norm = np.linalg.norm

#-----------------------------------
list_nx = [] ; list_ny = []
#list_nx += [1] ; list_ny += [1]
list_nx += [3] ; list_ny += [3]
list_nx += [7] ; list_ny += [7]
list_nx += [15] ; list_ny += [15]
#list_nx += [31] ; list_ny += [31]
#list_nx += [63] ; list_ny += [63]
#list_nx += [127] ; list_ny += [127]
#list_nx += [255] ; list_ny += [255]
#list_nx += [511] ; list_ny += [511]

nlevels = len(list_nx)

px      =  3  ; py      = 3

# ...
nnl         = 1
gamma       = 1
tol_ini     = 1.e-10
#tol_ini     = 5.e-7
#tol_ini     = 1.e-14
tol         = tol_ini
maxiter     = 1
maxiter_cg  = 6000

nu1         = 20
nu2         = 20
# ...
# ...
accel       = 'gmres'
#accel       = 'cg'
#maxiter_prec    = 5
maxiter_prec= 2
tol_prec    =1.e-10
# ...
# ...
maxiter_petsc  = 6000
tol_petsc  = 1.e-12
# ...
# ...
maxiter_pyamg = 200
tol_pyamg   = 1.e-10
# ...


# ...
#GalerkinCoarseGrid = True
GalerkinCoarseGrid = False
# ...

#withscipy   = True
withscipy   = False

withPETSc   = True
#withPETSc   = False

#IMPORT      = True
IMPORT      = False

#ALLSOLVERS  = True
ALLSOLVERS  = False
#-----------------------------------


#-----------------------------------
#           DIRICHLET
#-----------------------------------
# ...
#u_exact = lambda x,y : [exp ( 0.5 * ( x**2 + y**2 ) )]
#rho0 = lambda x,y : ( 1. + x**2 + y**2 ) * exp ( x**2 + y**2 )
#rho1 = lambda x,y : 1.
#c_rho = 1.
## ...
## ...
## values of u at the boundary
## ...
#AllDirichlet = True
# ...
#-----------------------------------

#-----------------------------------
#           MESH
#-----------------------------------
# ...
C0   = 1.0
rho0 = lambda x,y : 1.

C1   = 0.616805883732
t = 0.5
rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-t)**2+(y-0.5)**2-0.09)))

#    C1   =  1.75484181939
#    rho1 = lambda x,y : ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))

c_rho = C0/C1
# ...

# ...
# values of gradu.n at the boundary
# ...
def func_g(x,y):
    return [x,y]
# ...

# ...
bc_neumann={}
bc_neumann [0,0] = func_g
bc_neumann [0,1] = func_g
bc_neumann [0,2] = func_g
bc_neumann [0,3] = func_g
# ...

# ...
tc = {}
tc['A'] = lambda x,y : [1., 0., 0., 1.]
tc['b'] = lambda x,y : [1.e-3]
tc['u'] = lambda x,y : [0.]
tc['f'] = lambda x,y : [0.]
tc['bc_neumann'] = bc_neumann
# ...

AllDirichlet=True
#-----------------------------------

#-----------------------------------
list_geometry = []
for nx, ny in zip(list_nx, list_ny):
    geo = domain(n=[nx,ny], p=[px,py])
    list_geometry.append(geo)
#-----------------------------------

# ...

# ...
list_PDE = []
for geo in list_geometry:
#    PDE = picard(geometry=geo, bc_neumann=bc_neumann)
    PDE = picard(geometry=geo, AllDirichlet=AllDirichlet)
#    PDE = picard(geometry=geo, testcase=tc)
    PDE.meanConstraint = False
    list_PDE.append(PDE)
# ...

# ...
list_iPDE = []
if GalerkinCoarseGrid:
    list_iPDE = [-1]
else:
    list_iPDE = range(0, len(list_geometry))
# ...

# ...
geo_h = list_geometry[-1]
PDE_h = list_PDE[-1]

PDE_H = list_PDE[0]
# ...

# ...
def coarse_solver(b,x0=None):
    from scipy.sparse.linalg import cg
    A = PDE_H.system.get()
    return cg(A, b)[0]

#def coarse_solver(b,x0=None):
#    PDE_H.unknown.set(b)
#    PDE_H.solve(rho0, rho1, c_rho=c_rho, u0=x0 \
#                , maxiter=100, rtol=1.e-6, verbose=False)
#    return PDE_H.unknown.get()
# ...

# ...
MG = agregation(list_geometry, gamma, nu1, nu2 \
               , withscipy=withscipy \
               , withPETSc=withPETSc \
               , coarse_solver=coarse_solver)
# ...

# ...
print "Assembling "
for i in list_iPDE:
    PDE = list_PDE[i]
    i_ = i
    if i_<0:
        i_ = len(list_iPDE) + i
    PDE.assembly()
    PDE.system.save("A_"+str(i)+".mtx")
print "done."
# ...

# ...
list_A = []
for i in list_iPDE:
    PDE = list_PDE[i]
    list_A.append(PDE.system.get())
# ...

# ...
print "initializing MG "
MG.initialize(list_A)
print "done."
# ...

print "---------------------"
print "nx, ny   : ", list_nx, list_ny
print "px, py   : ", px, py
print "nlevels  : ", nlevels
print "gamma    : ", gamma
print "tol      : ", tol_ini
print "maxiter_cg : ", maxiter_cg
print "nu1 : ", nu1
print "nu2 : ", nu2
print "---------------------"



# ----------------------------------------------
def iterMG(b, x0, tol, maxiter):
    x = np.zeros_like(x0)
    mg_residuals = []

    t_start = time.time()
    x += MG.solve(b, x0=x0, tol=tol, maxiter=maxiter, residuals=mg_residuals)
    t_end = time.time()
    mg_elapsed = t_end - t_start

    # Compute (geometric) convergence factors
    mg_factor = (mg_residuals[-1]/mg_residuals[0])**(1.0/len(mg_residuals))
    # Compute relative residuals
    mg_residuals= np.array(mg_residuals)/mg_residuals[0]

    PDE_h.unknown.set(x)

#    errNorm = PDE_h.norm(exact=u_exact)

    etiq = '_p' + str(px) + '_n' + str(nx)
    np.savetxt('residuals'+etiq+'.txt', np.array(mg_residuals))
    np.savetxt('residuals.txt', np.array(mg_residuals))
    print "mg elapsed time  : ", mg_elapsed
    print "mg converges after ", MG.nloop, " V-cycles with error : ",mg_residuals[-1]
    print "mg convergence factor: %g"%(mg_factor)
    print "final error : ", np.linalg.norm(b-MG.A.dot(x))
#    print "exact error : ", errNorm

    return x
# ----------------------------------------------

# ----------------------------------------------
#       Compute an initial guess
# ----------------------------------------------
PDE_H.solve(rho0, rho1, c_rho=c_rho, u0=None \
            , maxiter=100, rtol=1.e-6, verbose=False)
vH = PDE_H.unknown.get()
vh = vH
for ilvl in range(0, nlevels):
    vh = MG.interpolation(ilvl, vh)
x0 = vh
# ----------------------------------------------

# ----------------------------------------------
print "---------------------"
print "MG-standalone solve"
# ...
U_h = PDE_h.unknown
rhs_h = PDE_h.rhs

U_h.set(vh)
PDE_h.update()

b = rhs_h.get()
x = np.zeros_like(b)

for i in range(0, nnl):
    print ">>"
    x = iterMG(b, x0, tol, maxiter)
    PDE_h.unknown.set(x)
    PDE_h.update()
    b  = PDE_h.rhs.get()
    print ">>>> error ", norm(x-x0)
    print x0[:5]
    print x[:5]
    x0 = x
print "done."
print "---------------------"
# ...
# ----------------------------------------------

## ...
#print ">>> Solving using Picard <<<"
## ...
#if PDE.Dirichlet:
#    U = PDE.unknown_dirichlet
#else:
#    U = PDE.unknown
## ...
#
## ...
#C0   = 1.0
#rho0 = lambda x,y : 1.
#
#C1   = 0.616805883732
#t = 0.5
#rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-t)**2+(y-0.5)**2-0.09)))
#
##    C1   =  1.75484181939
##    rho1 = lambda x,y : ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
#
#c_rho = C0/C1
## ...
#
## ...
#V = PDE.space
#V.nderiv_pts = 2
## ...
#
## ...
#PDE.solve(rho0, rho1, c_rho=None, u0=None, maxiter=100, rtol=1.e-6, verbose=True)
## ...
#
## ...
#PDE.plotMesh(ntx=60, nty=60)
## ...

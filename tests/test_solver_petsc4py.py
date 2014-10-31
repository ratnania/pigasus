# -*- coding: UTF-8 -*-
#! /usr/bin/python
from pigasus.utils.manager import context

try:
    from petsc4py import PETSc
    with_PETSC = True
except ImportError:
    with_PETSC = False
import numpy as np
import time
from scipy import sparse
from pigasus.gallery.poisson import poisson
from caid.cad_geometry import square as domain
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

tol_petsc = 1.e-12
maxiter_petsc = 6000

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

if with_PETSC:
    with context():

        PDE = poisson(geometry=geo, bc_dirichlet=bc_dirichlet, bc_neumann=bc_neumann,
                      AllDirichlet=AllDirichlet, metric=Metric)
        PDE.assembly()
        PDE.solve()

        # getting scipy matrix
        A_scipy = PDE.system.get()
        PDE.system.save("sys.mtx")

        b = np.ones(PDE.size)

        A = PETSc.Mat().createAIJ(size=A_scipy.shape,csr=(A_scipy.indptr, A_scipy.indices, A_scipy.data))
        # ...

        # Initialize ksp solver.
        ksp = PETSc.KSP().create()
        ksp.setOperators(A)
        pc = ksp.getPC()

        # Allow for solver choice to be set from command line with -ksp_type <solver>.
        # Recommended option: -ksp_type preonly -pc_type lu
        ksp.setFromOptions()

        ksptype = PETSc.KSP.Type.CG             ; pctype = None
        #ksptype = PETSc.KSP.Type.GMRES          ; pctype = None
        #ksptype = PETSc.KSP.Type.BICG           ; pctype = None
        #ksptype = PETSc.KSP.Type.BCGS           ; pctype = None
        #ksptype = PETSc.KSP.Type.BCGSL          ; pctype = None
        #ksptype = PETSc.KSP.Type.CGS            ; pctype = None
        #ksptype = PETSc.KSP.Type.STCG           ; pctype = None
        #ksptype = PETSc.KSP.Type.CGNE           ; pctype = None
        #ksptype = PETSc.KSP.Type.RICHARDSON     ; pctype = None
        #ksptype = PETSc.KSP.Type.CHEBYSHEV      ; pctype = None
        #ksptype = PETSc.KSP.Type.FGMRES         ; pctype = None
        #ksptype = PETSc.KSP.Type.CG             ; pctype = pc.Type.GAMG
        #ksptype = PETSc.KSP.Type.GMRES          ; pctype = pc.Type.GAMG

        ksp.setType(ksptype)
        if pctype is not None:
            pc.setType(pctype)

        txt = 'Using Petsc-' + str(ksp.getType())
        if pctype is not None:
            txt += ' with ' + str(pctype) + ' preconditioner.'
        print(txt)

        n,m = A.getSize()
        x0 = np.zeros(n)

        _b = PETSc.Vec().createWithArray(b, comm=PETSc.COMM_SELF)
        _x = PETSc.Vec().createWithArray(x0, comm=PETSc.COMM_SELF)

        ksp.rtol = tol_petsc
        ksp.max_it = maxiter_petsc
        ksp.setConvergenceHistory()

        # Solve!
        t_start = time.time()
        try:
            ksp.solve(_b, _x)
        except:
            print("PETCs-" + str(ksp.getType())+" Couldn't converge")

        t_end = time.time()
        elapsed = t_end - t_start

        r = _b.duplicate()
        u = _x.duplicate()
        ksp.buildSolution(u)
        ksp.buildResidual(u)
        try:
            err = ksp.getConvergenceHistory()[-1]
        except:
            err = None

        petsc_elapsed = 'Petsc-'+ str(ksp.getType())
        if pctype is not None:
            petsc_elapsed += ' with ' + str(pctype) + ' PC'
        petsc_elapsed += " \t : " + str(elapsed)

        petsc_txt = 'Petsc '+ str(ksp.getType())
        if pctype is not None:
            petsc_txt += ' with ' + str(pctype) + ' PC '
        petsc_txt += ' Converges after ' + str(ksp.getIterationNumber()) + ' iterations with error : ' + str(err)
        x = _x.getArray()
        petsc_txt += "\n final error : "+ str(np.linalg.norm(b-A_scipy * x))

        if ksptype == PETSc.KSP.Type.CG:
            petsc_cg_tol = err

        print("Elapsed time ", petsc_elapsed)
        print(petsc_txt)

        PDE.free()

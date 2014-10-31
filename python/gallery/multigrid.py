# -*- coding: UTF-8 -*-
from pigasus.multigrid.agregation import *
import numpy                as np
import sys
import time

__all__ = ['multigrid']

class multigrid(object):
    def __init__(self, PDE, list_geometry=None \
                 , gamma=1, nu1=10, nu2=10 \
                 , GalerkinCoarseGrid=True \
                 , withscipy=False, withPETSc=False):
        """
        list_geometry does not contain the finest geometry
        """
        self.PDE_h = PDE
        self.geo_h = PDE.geometry
        self.testcase = PDE.testcase

        self.nu1 = nu1
        self.nu2 = nu2
        self.gamma = gamma
        self.GalerkinCoarseGrid = GalerkinCoarseGrid
        self.withscipy = withscipy
        self.withPETSc = withPETSc

        self.initialized = False

        if list_geometry is None:
            print("You must give list_geometry.")
            print("Pigasus will stop.")
            sys.exit(0)

        # ...
        MG = agregation(list_geometry, gamma, nu1, nu2 \
                       , withscipy=withscipy, withPETSc=withPETSc)
        self.MG = MG
        # ...

    def initialize(self, list_A=None, verbose=False):
        # ...
        if self.initialized:
            return
        else:
            self.initialized = True
        # ...

        # ...
        list_DirFaces = self.PDE_h.list_DirFaces
        GalerkinCoarseGrid = self.GalerkinCoarseGrid
        # ...

        # ...
        if verbose:
            print("initializing MG ")
        A_h = self.PDE_h.system
        self.MG.initialize(A_h, list_A=list_A \
              , list_DirFaces=list_DirFaces \
              , GalerkinCoarseGrid=GalerkinCoarseGrid)
        if verbose:
            print("done.")
        # ...

    def solve(self, rhs, verbose=False, maxiter=200, tol=1.e-10 \
              , accel=None, maxiter_prec=3, tol_prec=1.e-10):

        self.initialize()

        # ----------------------------------------------
        # ...
        b = rhs.get()
        x = np.zeros_like(b)
        mg_residuals = []

        if accel is None:
            t_start = time.time()
            x += self.MG.solve(b, tol=tol, maxiter=maxiter, residuals=mg_residuals)
            t_end = time.time()
            mg_elapsed = t_end - t_start

            # Compute (geometric) convergence factors
            mg_factor = (mg_residuals[-1]/mg_residuals[0])**(1.0/len(mg_residuals))
            # Compute relative residuals
            mg_residuals= np.array(mg_residuals)/mg_residuals[0]

            if verbose:
                print("=======================================")
                print("==                                   ==")
                print("==       MG-standalone solve         ==")
                print("==                                   ==")
                print("=======================================")
                print("mg elapsed time  : ", mg_elapsed)
                print("mg converges after ", self.MG.nloop, " V-cycles with error : ",mg_residuals[-1])
                print("mg convergence factor: %g"%(mg_factor))
                print("final error : ", np.linalg.norm(b-self.MG.A.dot(x)))
                print("=======================================")
        else:
            t_start = time.time()
            x = self.MG.solve(b, tol=tol, maxiter=maxiter \
                         , maxiter_prec=maxiter_prec, tol_prec=tol_prec \
                         , accel=accel, residuals=mg_residuals)
            t_end = time.time()

            mgaccel_elapsed = t_end - t_start
            mgaccel_err = mg_residuals[-1]
            # Compute (geometric) convergence factors
            mgaccel_factor = (mg_residuals[-1]/mg_residuals[0])**(1.0/len(mg_residuals))
            # Compute relative residuals
            mg_residuals= np.array(mg_residuals)/mg_residuals[0]

            if verbose:
                print("=======================================")
                print("==                                   ==")
                print("==       MG-accelerated solve        ==")
                print("==                                   ==")
                print("=======================================")
                print("mgaccel elapsed time  : ", mgaccel_elapsed)
                print("mgaccel converges after ", len(mg_residuals), " iterations with error : ", mg_residuals[-1])
                print("mgaccel convergence factor: %g"%(mgaccel_factor))
                print("final error : ", np.linalg.norm(b-self.MG.A.dot(x)))
                print("=======================================")

        return mg_residuals
        # ...


if __name__ == '__main__':
    from caid.cad_geometry import square as domain
    from pigasus.gallery.poisson import *

    px = 2 ; py = 2

    list_nx = [] ; list_ny = []
    list_nx += [3] ; list_ny += [3]
    list_nx += [7] ; list_ny += [7]
    list_nx += [15] ; list_ny += [15]
    list_nx += [31] ; list_ny += [31]

    list_geometry = []
    for (nx,ny) in zip(list_nx, list_ny):
        geo = domain(n=[nx,ny], p=[px,py])
        list_geometry.append(geo)

    geo_h = list_geometry[-1]
    PDE = poisson(geometry=geo_h, AllDirichlet=True)
    rhs = PDE.rhs

    MG = multigrid(PDE, list_geometry=list_geometry)
    PDE.assembly()

    mg_residuals = MG.solve(rhs, verbose=True, accel=None)
#    mg_residuals = MG.solve(rhs, verbose=True, accel='gmres')

# -*- coding: UTF-8 -*-
#! /usr/bin/python

from caid.cad_geometry import square
from caid.cad_geometry import circle
from caid.cad_geometry import quart_circle
from caid.cad_geometry import annulus
from matplotlib import pyplot as plt
import numpy                as np
from time import time
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)

# ... import picard from monge_ampere module
from pigasus.utils.load import load
monge_ampere    = load("monge_ampere")
picard          = monge_ampere.picard
# ...

abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)


class adaptiveMeshMA(object):
    def __init__(self, geo_h, geo_H=None, verbose=False,
                 bc_neumann=None):
        # ...
        if geo_H is None:
            withTwoGrids = False
        else:
            withTwoGrids = True
        # ...

        if withTwoGrids:
            PDE_H = picard(geometry=geo_H, bc_neumann=bc_neumann)
        else:
            PDE_H = None
        PDE_h = picard(geometry=geo_h, bc_neumann=bc_neumann)

        self.withTwoGrids = withTwoGrids
        self.verbose = verbose
        self.geo_H = geo_H
        self.geo_h = geo_h
        self.PDE_H = PDE_H
        self.PDE_h = PDE_h

    def construct(self, rho0, rho1, c_rho=None \
                  , rtol_h = 1.e-3, rtol2_h = 1.e-9, maxiter_h = 40 \
                  , rtol_H = 1.e-5, rtol2_H = 1.e-6, maxiter_H = 40 \
                 ):
        withTwoGrids = self.withTwoGrids
        verbose = self.verbose
        PDE_h = self.PDE_h
        PDE_H = self.PDE_H
        geo_h = self.geo_h
        geo_H = self.geo_H

        # ...
        if withTwoGrids:
            if PDE_H.Dirichlet:
                U_H = PDE_H.unknown_dirichlet
            else:
                U_H = PDE_H.unknown

        if PDE_h.Dirichlet:
            U_h = PDE_h.unknown_dirichlet
        else:
            U_h = PDE_h.unknown
        # ...

        # ...
        if withTwoGrids:
            if verbose:
                print("*****************************")
                tb = time()

            Errors_H, ErrorsH1_H = PDE_H.solve(  rho0, rho1, c_rho=c_rho, u0=None \
                        , maxiter=maxiter_H, rtol=rtol_H, rtol2=rtol2_h, verbose=verbose)

            if verbose:
                te = time()
                print("Coarse solver converges after ", len(Errors_H) \
                        , " with final error ", Errors_H[-1] \
                        , " with final H1-error ", ErrorsH1_H[-1])
                print("Elapsed time ", te-tb)
                print("*****************************")

            PDE_H.transferSolution(geo_H, U_H, geo_h, U_h)
            u0 = U_h.get()
        else:
            u0 = np.zeros(PDE_h.size)
        # ...

        # ...
        if verbose:
            print("*****************************")
            tb = time()

        Errors_h, ErrorsH1_h = PDE_h.solve(  rho0, rho1, c_rho=None, u0=u0 \
                    , maxiter=maxiter_h, rtol=rtol_h, rtol2=rtol2_h, verbose=verbose)

        if verbose:
            te = time()
            print("Monge-Ampere eq. converges after ", len(Errors_h) \
                    , " with final error ", Errors_h[-1] \
                    , " with final H1-error ", ErrorsH1_h[-1])
            print("Elapsed time ", te-tb)
            print("*****************************")
        # ...

        # ...
        if withTwoGrids:
            uH = U_H.get()
        uh = U_h.get()

        if verbose:
            if withTwoGrids:
                print("Error-coarse        ", np.abs(1.-PDE_H.norm(exact=PDE_H.Err_func)))
            print("Error-fine          ", np.abs(1.-PDE_h.norm(exact=PDE_h.Err_func)))

        if withTwoGrids:
            U_H.set(uH)
        U_h.set(uh)
        # ...

        # ...
        return self.genGeometry()
        # ...

    def plot(self):
        PDE_h = self.PDE_h
        PDE_h.plotMesh(ntx=60, nty=60)

    def genGeometry(self):
        from caid.op_nurbs import opNURBS, grad
        from caid.cad_geometry import cad_geometry, cad_grad_nurbs, cad_nurbs

        geo = self.geo_h
        PDE = self.PDE_h
        if PDE.Dirichlet:
            U = PDE.unknown_dirichlet
        else:
            U = PDE.unknown
        geo_f = cad_geometry()
        for patch_id in range(0, geo.npatchs):
            met   = geo[patch_id]

            _C = U.tomatrix(patch_id)
            C = np.zeros_like(met.points)
            C[...,:2] = met.points[...,:2]
            C[...,2] = _C
            srf = cad_nurbs(met.knots, C, weights= met.weights)
            gsrf = grad(srf)
            grad_nrb = cad_grad_nurbs(gsrf)
            print("name ", grad_nrb.__class__.__name__)
            geo_f.append(grad_nrb)
            geo_f[-1].set_orientation(met.orientation)
            geo_f[-1].set_rational(met.rational)

        geo_f.set_internal_faces(geo.internal_faces)
        geo_f.set_external_faces(geo.external_faces)
        geo_f.set_connectivity(geo.connectivity)

        return geo_f

#-----------------------------------


#-----------------------------------
if __name__ == '__main__':
    verbose = False
    withTwoGrids = True

    #p_H = [ 5, 5]
    #p_h = [ 5, 5]

    p_H = [ 3, 3]
    p_h = [ 3, 3]

    #p_H = [ 2, 2]
    #p_h = [ 2, 2]

    rtol_H = 1.e-3
    rtol_h = 1.e-3

    rtol2_H = 1.e-6
    rtol2_h = 1.e-9

    maxiter_H = 40
    maxiter_h = 40

    n_H = [3,3]
    n_h = [7,7]
    #n_h = [15,15]
    #-----------------------------------

    #-----------------------------------
    # ...
    # values of gradu.n at the boundary
    # ...
    def func_g(x,y):
        return [x,y]
    # ...

    # ...
    # values of u at the boundary
    # ...
    bc_neumann={}
    bc_neumann [0,0] = func_g
    bc_neumann [0,1] = func_g
    bc_neumann [0,2] = func_g
    bc_neumann [0,3] = func_g
    # ...

    # ...
    # exact solution
    # ...
    C0 = 1.0

    # ... test 1
    rho0 = lambda x,y : 1.
    C1   = 0.616805883732
    t = 0.5
    rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-t)**2+(y-0.5)**2-0.09)))
    # ...

    # ...
    if withTwoGrids:
#        geo_H = square(n=n_H, p=p_H)
        geo_H = circle(radius=2.0, n=n_H, p=p_H)
        P = geo_H[0].points
        P[...,:2] += 1.
        geo_H[0].set_points(P)
#    geo_h = square(n=n_h, p=p_h)
    geo_h = circle(radius=2.0, n=n_h, p=p_h)
    P = geo_h[0].points
    P[...,:2] += 1.
    geo_h[0].set_points(P)
    # ...

    # ...
    adaptMesh = adaptiveMeshMA(geo_h, geo_H=geo_H, verbose=verbose, bc_neumann=bc_neumann)
    # ...

    # ...
    geo_f = adaptMesh.construct(rho0, rho1)
    geo_f.save("adaptive_mesh.xml")
    adaptMesh.plot()
    plt.show()
    # ...

    from pigasus.fem.metric import metric
    Metric = metric(geometry=geo_f)



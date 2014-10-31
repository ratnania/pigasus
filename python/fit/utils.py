# -*- coding: UTF-8 -*-
#! /usr/bin/python

import matplotlib.pyplot    as plt
import numpy                as np
from time import time
from pigasus.fem.basicPDE import basicPDE
from caid.cad_geometry import cad_geometry, cad_nurbs
from caid.cad_geometry import square as patch
from caid.core.bspline import bsp
#from igakit.nurbs import NURBS
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

#-----------------------------------
# ...
def compute_uk(list_Q, method="chord"):
    """
    this routine computes the parameters uk st C(uk) = Vk
    """
    npts = len(list_Q)
    lpr_u = np.zeros(npts)
    lpr_u[0]  = 0.
    lpr_u[-1] = 1.

    def compute_distance(A,B):
        return np.sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)

    if method == "uniform":
        lpr_u = np.linspace(0,1,npts)
        return lpr_u

    if method == "chord":
        d = 0.
        for (Qi,Qi1) in zip(list_Q[0:-1], list_Q[1:]):
            d += compute_distance(Qi, Qi1)
        for k in range(1,npts-1):
            Qi = list_Q[k] ; Qi1 = list_Q[k+1]
            dk = compute_distance(Qi, Qi1)
            lpr_u[k] = lpr_u[k-1] + dk / d
        return lpr_u

    if method == "centripetal":
        d = 0.
        for (Qi,Qi1) in zip(list_Q[0:-1], list_Q[1:]):
            d += np.sqrt(compute_distance(Qi, Qi1))
        for k in range(1,npts-1):
            Qi = list_Q[k] ; Qi1 = list_Q[k+1]
            dk = compute_distance(Qi, Qi1)
            lpr_u[k] = lpr_u[k-1] + np.sqrt(dk) / d
        return lpr_u
# ...


# ...
def evalBasis1D(nrb, list_uk, rational=0):
    Ux = nrb.knots[0]
    nx, = nrb.shape
    px, = nrb.degree
    nen = (px+1)
    Ww = np.array(nrb.weights)
    list_basis = []
    for uk in list_uk:
        C = bsp.AssembleBasis1(0,0, rational,px,Ux,Ww,np.array([uk]))
        list_basis.append(C[0,:,0]) # le premier indice de C est la taille de uk
    return list_basis
# ...

# ...
def addContributions1D(nrb, list_uk, list_basis, ID, Mat_shape):
    """
    """
    rows = [] ; cols = [] ; data = []
    Ux = nrb.knots[0]
    nx, = nrb.shape
    px, = nrb.degree
    nen = (px+1)
    for ind in range(0,len(list_uk)):
        uk    = list_uk[ind]
        basis = list_basis[ind]
        spanx = bsp.FindSpan(px,Ux,uk)
        # b = by * (px+1) + bx
        for b1, v1 in enumerate(basis):
            bx = b1
            Ix = spanx - px + bx
            A1 = ID[Ix]

            if (A1 > 0):
                rows.append(A1-1)
                cols.append(ind)
                data.append(v1)

    rows = np.array(rows)
    cols = np.array(cols)
    data = np.array(data)

    shp = [Mat_shape[0],len(list_uk)]

    D   = coo_matrix((data, (rows, cols)), shape=shp)
    D   = D.tocsr()
    Dt  = D.transpose().tocsr()
    DtD = D * Dt
    return D, DtD.transpose().tocsr()
# ...

# ...
def computeC1Coef1D(nrb,face):
    v, = nrb.knots
    q, = nrb.degree
    if face == 0:
        k = q+1
        d = v[k]-v[0]
    if face == 1:
        k = q+1
        d = v[-1] - v[-(k+1)]

    c = (k-1) / d
    return c
# ...

# ...
def genLineC0Constraint(  Rd, list_Am \
                          , list_As=None \
                          , list_vals=None \
                         ):
    rows = [] ; values = [] ; rhs = []
    if list_As is not None:
        for (Am, As) in zip(list_Am, list_As):
            rows.append([Am, As])
            values.append([1., -1.])
            rhs.append([0.]*Rd)
    if list_vals is not None:
        for (Am, vals) in zip(list_Am, list_vals):
            rows.append([Am])
            values.append([1.])
            rhs.append(vals)

    return rows, values, rhs
# ...

# ...
def genLineC1Constraint(  Rd, list_Am, list_Bm, cm \
                          , list_As=None, list_Bs=None, cs=None \
                          , list_vals=None \
                         ):
    rows = [] ; values = [] ; rhs = []
    if (list_As is not None) and (list_Bs is not None):
        for (Am, Bm, As, Bs) in zip(list_Am, list_Bm, list_As, list_Bs):
            rows.append([Am, Bm, As, Bs])
            values.append([-cm, cm, -cs, cs])
            rhs.append([0.]*Rd)
    if (list_vals is not None) :
        for (Am, Bm, vals) in zip(list_Am, list_Bm, list_vals):
            rows.append([Am, Bm])
            values.append([-cm, cm])
            rhs.append(vals)

    return rows, values, rhs
# ...
#-----------------------------------
# ...
def evalBasis2D(nrb, list_uk, list_vk \
                , rational=0, verbose=False):
    if verbose:
        tb = time()
    Ux = nrb.knots[0]
    Uy = nrb.knots[1]
    nx, ny = nrb.shape
    px, py = nrb.degree
    nen = (px+1)*(py+1)
    Ww = np.array(nrb.weights)
    list_basis = []
    for (uk,vk) in zip(list_uk, list_vk):
        C = bsp.AssembleBasis2(0,0,rational,px,Ux,py,Uy,Ww,np.array([uk]),np.array([vk]))
        list_basis.append(C[0,:,0]) # le premier indice de C est la taille de uk
    if verbose:
        te = time()
        print("elapsed time basis evaluation ", te-tb)
    return list_basis
# ...

# ...
def addContributions2D(nrb, list_uk, list_vk \
                       , list_basis, ID \
                       , verbose=False):
    """
    """
    if verbose:
        tb = time()
    rows = [] ; cols = [] ; data = []
    Ux = nrb.knots[0]
    Uy = nrb.knots[1]
    nx, ny = nrb.shape
    px, py = nrb.degree
    nen = (px+1)*(py+1)
    for ind in range(0,len(list_uk)):
        uk    = list_uk[ind]
        vk    = list_vk[ind]
        basis = list_basis[ind]
        spanx = bsp.FindSpan(px,Ux,uk)
        spany = bsp.FindSpan(py,Uy,vk)
        # b = by * (px+1) + bx
        for b1, v1 in enumerate(basis):
            by = b1 / (px+1)
            bx = b1 - by * (px+1)
            Ix = spanx - px + bx
            Iy = spany - py + by
            A1 = ID[Ix, Iy]

            if (A1 > 0):
                rows.append(A1-1)
                cols.append(ind)
                data.append(v1)

    rows = np.array(rows)
    cols = np.array(cols)
    data = np.array(data)

    if verbose:
        te = time()
        print("elapsed time discrete Mass ", te-tb)
    return rows, cols, data
# ...

# ...
def computeC1Coef2D(nrb,face):
    u,v = nrb.knots
    p,q = nrb.degree
    if face == 0:
        k = q+1
        d = v[k]-v[0]
    if face == 1:
        k = p+1
        d = u[k]-u[0]
    if face == 2:
        k = q+1
        d = v[-1] - v[-(k+1)]
    if face == 3:
        k = p+1
        d = u[-1] - u[-(k+1)]

    c = (k-1) / d
    return c
# ...

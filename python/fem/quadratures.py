# -*- coding: UTF-8 -*-
#! /usr/bin/python

"""
This module contains some routines to generate quadrature points in 1D
it has also a routine uniform, which generates uniform points
with weights equal to 1
"""

__author__="ARA"
__all__ = ['quadratures']
__date__ ="$Jan 13, 2012 11:55:18 AM$"

import numpy as np
class quadratures:
    def __init__(self):
        pass

    def uniform(self,k):
        xg = np.linspace(-1., 1., k+1)
        w  = np.ones(k+1)

        return xg, w

    def gauss_lobatto(self,k):
        beta = .5 / np.sqrt(1-(2 * np.arange(1., k + 1)) ** (-2)) #3-term recurrence coeffs
        beta[-1] = np.sqrt((k / (2 * k-1.)))
        T = np.diag(beta, 1) + np.diag(beta, -1) # jacobi matrix
        D, V = np.linalg.eig(T) # eigenvalue decomposition
        xg = np.real(D); i = xg.argsort(); xg.sort() # nodes (= Legendres points)
        w = 2 * (V[0, :]) ** 2; # weights

        return xg, w[i]

    def gauss_legendre(self,ordergl,tol=10e-9):
        ''' x,A = gaussNodes(m,tol=10e-9)
            Returns nodal abscissas {x} and weights {A} of
            Gauss-Legendre m-point quadrature.
        '''
        m = ordergl + 1
        from math import cos,pi
        from numpy import zeros

        def legendre(t,m):
            p0 = 1.0; p1 = t
            for k in range(1,m):
                p = ((2.0*k + 1.0)*t*p1 - k*p0)/(1.0 + k )
                p0 = p1; p1 = p
            dp = m*(p0 - t*p1)/(1.0 - t**2)
            return p1,dp

        A = zeros(m)
        x = zeros(m)
        nRoots = (m + 1)/2          # Number of non-neg. roots
        for i in range(nRoots):
            t = cos(pi*(i + 0.75)/(m + 0.5))  # Approx. root
            for j in range(30):
                p,dp = legendre(t,m)          # Newton-Raphson
                dt = -p/dp; t = t + dt        # method
                if abs(dt) < tol:
                    x[i] = t; x[m-i-1] = -t
                    A[i] = 2.0/(1.0 - t**2)/(dp**2) # Eq.(6.25)
                    A[m-i-1] = A[i]
                    break
        return x,A

    def generate(self,apr_a, k,as_type="lobatto"):
#    def generate(self,a, b, N, k,as_type="lobatto"):
        """
        this routine generates a quad pts on the grid linspace(a,b,N)
        """
        if k == 1:
            x = np.asarray([-1., 1.])
            w = np.asarray([0.5, 0.5])
            grid = apr_a
            N = len(apr_a)
            xgl = np.zeros((N-1, k + 1))
            wgl = np.zeros((N-1, k + 1))
            for i in range (0, N-1):
                xmin = grid[i];xmax = grid[i + 1];dx = 0.5 * (xmax-xmin)
                tab = dx * x + dx + xmin
                xgl[i, :] = tab[::-1]
                wgl[i, :] = 0.5 * ( xmax - xmin ) * w

            return xgl,wgl

        if as_type.split('-')[0] == "radau":
            # this will generates radau qd points at the left and right of the interval
            # second rule is needed for the inside
            grid = apr_a
            N = len(apr_a)
            xgl = np.zeros((N-1, k + 1))
            wgl = np.zeros((N-1, k + 1))

            # ...
            # left
            # ...
            import qd_radau as qr
            x, w = qr.radau_left(k)
            x = x[::-1] ; w = w[::-1]
            i = 0
            xmin = grid[i];xmax = grid[i + 1];dx = 0.5 * (xmax-xmin)
            tab = dx * x + dx + xmin
            xgl[i, :] = tab[::-1]
            wgl[i, :] = 0.5 * ( xmax - xmin ) * w
            # ...

            # ...
            # inside
            # ...
            if as_type.split('-')[1] == "uniform":
                x, w = self.uniform(k)
            if as_type.split('-')[1] == "lobatto":
                x, w = self.gauss_lobatto(k)
            if as_type.split('-')[1] == "legendre":
                x, w = self.gauss_legendre(k)
            for i in range (1, N-2):
                xmin = grid[i];xmax = grid[i + 1];dx = 0.5 * (xmax-xmin)
                tab = dx * x + dx + xmin
                xgl[i, :] = tab[::-1]
                wgl[i, :] = 0.5 * ( xmax - xmin ) * w
            # ...

            # ...
            # right
            # ...
            import qd_radau as qr
            x, w = qr.radau_right(k)
            x = x[::-1] ; w = w[::-1]
            i = N-2
            xmin = grid[i];xmax = grid[i + 1];dx = 0.5 * (xmax-xmin)
            tab = dx * x + dx + xmin
            xgl[i, :] = tab[::-1]
            wgl[i, :] = 0.5 * ( xmax - xmin ) * w
            # ...

            return xgl,wgl
        if as_type == "uniform":
            x, w = self.uniform(k)
        if as_type == "lobatto":
            x, w = self.gauss_lobatto(k)
        if as_type == "legendre":
            x, w = self.gauss_legendre(k)
        if as_type == "radau_left":
            import qd_radau as qr
            x, w = qr.radau_left(k)
            x = x[::-1] ; w = w[::-1]
        if as_type == "radau_right":
            import qd_radau as qr
            x, w = qr.radau_right(k)
            x = x[::-1] ; w = w[::-1]
#        grid = np.linspace(a, b, N)
        grid = apr_a
        N = len(apr_a)
        xgl = np.zeros((N-1, k + 1))
        wgl = np.zeros((N-1, k + 1))
        for i in range (0, N-1):
            xmin = grid[i];xmax = grid[i + 1];dx = 0.5 * (xmax-xmin)
            tab = dx * x + dx + xmin
            xgl[i, :] = tab[::-1]
            wgl[i, :] = 0.5 * ( xmax - xmin ) * w

        return xgl,wgl

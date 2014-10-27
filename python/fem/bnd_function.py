# -*- coding: UTF-8 -*-
import numpy as np
from scipy.interpolate import splev, splrep

# ...
class bnd_function(object):
    # ...
    def __init__(self, nrb, face, g, sites=None):
        """
        sites : the parametric 1D points where to compute g
        """
        self.nrb    = nrb #.clone().transpose([1,0])
#        print nrb.points
#        print self.nrb.points
        self.face   = face
        self.dim    = nrb.dim
        self.x_bnd  = None
        self.xi_bnd = None
        self.sites  = sites
        self.g      = g
        self.gi     = None
        self.axis   = None
        self.ind    = None

        # ...
        if self.dim == 1:
            self.axis = 0

            if face in [0]:
                self.ind = 0
            if face in [1]:
                self.ind = -1

            self.set_bnd_data_1D()
        # ...

        # ...
        if self.dim == 2:
            if face in [0,2]:
                self.axis = 0
            if face in [1,3]:
                self.axis = 1

            if face in [0,1]:
                self.ind = 0
            if face in [2,3]:
                self.ind = -1

            self.set_bnd_data_2D()
        # ...
    # ...

    # ...
    def rise(self, type):
        if self.dim == 1:
            return self.rise_1D(type)

        if self.dim == 2:
            return self.rise_2D(type)
    # ...

    # ...
    def set_bnd_data_1D(self):
        """
        sets all data for the boundary
        """
        # ...
        n = self.nrb.shape[self.axis]
        p = self.nrb.degree[self.axis]
        knots = self.nrb.knots[self.axis]

        xb = knots[0] ; xe = knots[-1]

        if self.face in [0]:
            # we start by defining an interpolated function over this face
            # we compute a list of points from the face where we'll compute g
            xi_bnd = xb
            self.xi_bnd = [xb]

        # ...
        # ...
        if self.face in [1]:
            # we start by defining an interpolated function over this face
            # we compute a list of points from the face where we'll compute g
            xi_bnd = xe
            self.xi_bnd = [xe]
        # ...
        pts_bnd = self.nrb(u=xi_bnd)

        x_bnd = pts_bnd[0]
        self.x_bnd = [x_bnd]

        # compute the corresponding values of g
        g_bnd = self.g(x_bnd)
        self.gi = g_bnd

    # ...
    def set_bnd_data_2D(self):
        """
        sets all data for the boundary
        """
        # ...
        n = self.nrb.shape[self.axis]
        p = self.nrb.degree[self.axis]
        knots = self.nrb.knots[self.axis]

        xb = knots[0] ; xe = knots[-1]
        yb = self.nrb.knots[self.dim-self.axis-1][0] ; ye = self.nrb.knots[self.dim-self.axis-1][-1]
        t = knots[p+1:-p-1]

        if self.face in [0,2]:
            # we start by defining an interpolated function over this face
            # we compute a list of points from the face where we'll compute g
            if self.sites is None:
                xi_bnd  = np.linspace(xb,xe,2*(n+p+1))
                self.sites = xi_bnd
            else:
                xi_bnd = self.sites
            if self.face == 0:
                eta_bnd = np.asarray(yb)
                self.xi_bnd = [xi_bnd, yb * np.ones_like(xi_bnd)]
            if self.face == 2:
                eta_bnd = np.asarray(ye)
                self.xi_bnd = [xi_bnd, ye * np.ones_like(xi_bnd)]
        # ...
        # ...
        if self.face in [1,3]:
            # we start by defining an interpolated function over this face
            # we compute a list of points from the face where we'll compute g
            if self.sites is None:
                eta_bnd  = np.linspace(xb,xe,2*(n+p+1))
                self.sites = eta_bnd
            else:
                eta_bnd = self.sites
            if self.face == 1:
                xi_bnd = np.asarray(yb)
                self.eta_bnd = [yb * np.ones_like(eta_bnd), eta_bnd]
            if self.face == 3:
                xi_bnd = np.asarray(ye)
                self.eta_bnd = [ye * np.ones_like(eta_bnd), eta_bnd]
        # ...
        pts_bnd = self.nrb(u=xi_bnd, v=eta_bnd)
        x_bnd = pts_bnd[:,0]
        y_bnd = pts_bnd[:,1]
        self.x_bnd = [x_bnd, y_bnd]
#        print "=============="
#        print self.face
#        print x_bnd
#        print y_bnd

        # compute the corresponding values of g
        g_bnd = self.g(x_bnd, y_bnd)
#        print "g_bnd ", g_bnd
        # compute the interpolated function gi
        # ...
        # TODO to change for vectorial functions
#        print "---------"
#        print "face ", self.face
#        print 'xb ', xb, ', xe ', xe, ', p ', p
        self.tck = splrep(self.sites, g_bnd[0], xb=xb, xe=xe, k=p, t=t, s=0)
        # ...

        self.gi = lambda s : splev(s, self.tck, der=0)
#        print "-- TCK --"
#        print self.tck[0]
#        print self.tck[1]
#        print self.tck[2]
#        import pylab as pl
#        pl.plot(self.sites,g_bnd[0],'or')
#        pl.plot(self.sites,self.gi(self.sites))
#        pl.show()
    # ...

    # ...
    def get_value_1D(self, axis, xi):
        axis = 0
        return self.basisFct(axis,self.ind,xi)
    # ...

    # ...
    def get_value_2D(self, axis, xi):
        if self.face in [0,2]:
            if axis==0:
                return self.gi(xi)
            if axis==1:
                return self.basisFct(axis,self.ind,xi)
        if self.face in [1,3]:
            if axis==1:
                return self.gi(xi)
            if axis==0:
                return self.basisFct(axis,self.ind,xi)
    # ...

    # ...
    def rise_1D(self, type):
        """
        this routine will creates a 1D function from a given function "g" defined on
        the face "face"
        """

        if type == "function":
            # ...
            def gh(xi):
                return self.get_value_1D(0, xi)
            # ...

        if type == "coeff":
            list_n = self.nrb.shape
            gh = np.zeros(list_n)
            coef = self.gi[0]

            if self.face == 0:
                gh[0]  = coef
            if self.face == 1:
                gh[-1] = coef

            gh = gh

        return gh
    # ...

    # ...
    def rise_2D(self, type):
        """
        this routine will creates a 2D function from a given function "g" defined on
        the face "face"
        """

        if type == "function":
            self.gi = lambda xi : splev(xi, self.tck, der=0)

            # ...
            def gh(xi,eta):
                return self.get_value_2D(0, xi) * self.get_value_2D(1, eta)
            # ...

        if type == "coeff":
            list_n = self.nrb.shape
            gh = np.zeros(list_n)
            coef = self.tck[1][0:list_n[self.axis]]
#            print "%%%%%%%%%%%%%%%%%%%%%%%%%"
#            print "face   ", self.face
#            print "n      ", list_n
#            print "axis   ", self.axis
#            print "tck[1] ", self.tck[1]
#            print "coef   ", coef

            # TODO : to remove
#            coef[0] = 0.
#            coef[-1] = 0.
            #

            if self.face == 0:
                gh[:,0]  = coef
            if self.face == 1:
                gh[0,:]  = coef
            if self.face == 2:
                gh[:,-1] = coef
            if self.face == 3:
                gh[-1,:] = coef

            gh = gh.transpose()

        return gh
    # ...

    # ...
    def basisFct(self, axis, i, x):
        """evaluate b-spline starting at node i at x for the axis direction"""
        u = self.nrb.knots[axis]
        n = self.nrb.shape[axis]
        p = self.nrb.degree[axis]

        c = np.zeros_like(u)

        if i < 0:
            c[n + i]=1.
        else:
            c[i]=1.
        return splev(x,(u,c,p))

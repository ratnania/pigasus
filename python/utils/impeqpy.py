# -*- coding: UTF-8 -*-
#!/usr/bin/env python
from math import *
from numpy import *

"""
def f(x):
return x**2

def df(x):
return 2.0*x

def g(x,y):
return x**2+y**2-1

def dxg(x,y):
return 2.0*x

def dyg(x,y):
return 2.0*y
"""

class impeqpy:
    #this class contains some useful function for solving equation : f(x,y)
    def __init__ ( self, tol=1.e-9, maxniter = 5000, verbose=False ):
        self.tol = tol
        self.verbose = verbose
        self.maxniter = maxniter

    """
    Ubiquitous Newton-Raphson algorithm for solving
    f(x) = 0
    where a root is repeatedly estimated by
    x = x - f(x)/f'(x)
    until |dx|/(1+|x|) < TOL is achieved.  This termination condition is a
    compromise between
    |dx| < TOL,  if x is small
    |dx|/|x| < TOL,  if x is large
    """
    def newton(self, func, funcd, level,x=1.0):
        # f(x)=func(x), f'(x)=funcd(x)
        f, fd = func(x), funcd(x)
        f = f - level
        count = 0
        while (count < self.maxniter):
            dx = f / float(fd)
            if abs(dx) < self.tol * (1 + abs(x)):
                #print "newton(%d): x=%s, f(x)=%s" % (count, x, f)
                return x - dx
            x = x - dx
            f, fd = func(x), funcd(x)
            f = f - level
            count = count + 1
            if self.verbose and (count == self.maxniter):
                print("Warning: newton solver iterations = ", self.maxniter)
            #print "newton(%d): x=%s, f(x)=%s" % (count, x, f)

    def newton2D(self,func,dfunc,param,val, level=0.0,init=1.0):
        #dg may be dfunc/dx or dfunc/dy
        if (param==1):
            x=init
            y=val
        if (param==2):
            x=val
            y=init
        f, fd = func(x,y), dfunc(x,y)
        f = f - level
        count = 0
        #finding x
        if (param==1):
            while (count < self.maxniter):
                dx = f / float(fd)
                if abs(dx) < self.tol * (1 + abs(x)):
                    #print "newton(%d): x=%s, f(x)=%s" % (count, x, f)
                    return x - dx
                x = x - dx
                f, fd = func(x,y), dfunc(x,y)
                f = f - level
                count = count + 1
                if self.verbose and (count == self.maxniter):
                    print("Warning: newton solver iterations = ", self.maxniter)
                #print "newton(%d): x=%s, f(x)=%s" % (count, x, f)
                #if maxniter achieved
            return None
        #finding y
        if (param==2):
            while (count < self.maxniter):
                dy = f / float(fd)
                if abs(dy) < self.tol * (1 + abs(y)):
                    #print "newton(%d): x=%s, f(x)=%s" % (count, x, f)
                    return y - dy
                y = y - dy
                f, fd = func(x,y), dfunc(x,y)
                f = f - level
                count = count + 1
                if self.verbose and (count == self.maxniter):
                    print("Warning: newton solver iterations = ", self.maxniter)
                #print "newton(%d): x=%s, f(x)=%s" % (count, x, f)
                #if maxniter achieved
            return  None

    def solve2Dx(self,func,dfunc,level,apr_x,apr_y, y0=None):
        #apr_x input
        #apr_y output
        #apr_x array of size n
        n = size(apr_x)
        param = 2
        init = 1.
        for i in range(0,n):
            val = apr_x[i]
            #prendre l'initialisation par defaut
            if y0 is not None:
                init = y0[i]
            apr_y[i]=self.newton2D(func,dfunc,param,val, level, init=init)
            #specifier l'initialisation
            #init = apr_init[i]
            #self.newton2D(func,dfunc,param,val,init)

    def solve2Dy(self,func,dfunc,level,apr_x,apr_y, x0=None):
        #apr_x input
        #apr_y output
        #apr_x array of size n
        n = size(apr_y)
        param = 1
        init = 1.
        for i in range(0,n):
            val = apr_y[i]
            #prendre l'initialisation par defaut
            if x0 is not None:
                init = x0[i]
            apr_x[i]=self.newton2D(func,dfunc,param,val, level, init=init)
            #specifier l'initialisation
            #init = apr_init[i]
            #self.newton2D(func,dfunc,param,val,init)

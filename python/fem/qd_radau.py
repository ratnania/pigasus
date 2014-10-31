# -*- coding: UTF-8 -*-
import numpy as np
import scipy as sp
import scipy.linalg


def radau_left(n):
        if n == 0 :
            x = np.asarray([-1.0]) ; w = np.asarray([1.0])
            return x,w
	alpha,beta = __rec_jacobi__(n+1,0,1)
	x,w = __radau__(alpha,beta,-1)
        x[0] = -1.0
	return x,w

def radau_right(n):
        if n == 0 :
            x = np.asarray([1.0]) ; w = np.asarray([1.0])
            return x,w
	alpha,beta = __rec_jacobi__(n+1,0,1)
	x,w = __radau__(alpha,beta,1)
        x[-1] = 1.0
	return x,w

def __gauss__(alpha,beta):
    """
        Compute the Gauss nodes and weights from the recursion
        coefficients associated with a set of orthogonal polynomials

        Inputs:
        alpha - recursion coefficients
        beta - recursion coefficients

        Outputs:
        x - quadrature nodes
        w - quadrature weights

        Adapted from the MATLAB code by Walter Gautschi
        http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m
    """

    from scipy.linalg import eig_banded

    A = np.vstack((np.sqrt(beta),alpha))
    x,V = eig_banded(A,lower=False)
    w = beta[0]*sp.real(sp.power(V[0,:],2))
    return x,w

def __radau__(alpha,beta,xr):
    """
        Compute the Radau nodes and weights with the preassigned node xr

        Inputs:
        alpha - recursion coefficients
        beta - recursion coefficients
        xr - assigned node location

        Outputs:
        x - quadrature nodes
        w - quadrature weights

        Based on the section 7 of the paper "Some modified matrix eigenvalue
        problems" by Gene Golub, SIAM Review Vol 15, No. 2, April 1973, pp.318--334
    """
    from scipy.linalg import solve_banded
    n = len(alpha)-1
    f = np.zeros(n)
    f[-1] = beta[-1]
    A = np.vstack((np.sqrt(beta),alpha-xr))
    J = np.vstack((A[:,0:-1],A[0,1:]))
    delta = solve_banded((1,1),J,f)
    alphar = alpha
    alphar[-1] = xr+delta[-1]
    x,w = __gauss__(alphar,beta)
    return x,w

def __lobatto__(alpha,beta,xl1,xl2):
    """
        Compute the Lobatto nodes and weights with the preassigned node xl1,xl2

        Inputs:
        alpha - recursion coefficients
        beta - recursion coefficients
        xl1 - assigned node location
        xl2 - assigned node location

        Outputs:
        x - quadrature nodes
        w - quadrature weights

        Based on the section 7 of the paper "Some modified matrix eigenvalue
        problems" by Gene Golub, SIAM Review Vol 15, No. 2, April 1973, pp.318--334
    """
    from scipy.linalg import solve_banded, solve
    n = len(alpha)-1
    en = np.zeros(n)
    en[-1] = 1
    A1 = np.vstack((np.sqrt(beta),alpha-xl1))
    J1 = np.vstack((A1[:,0:-1],A1[0,1:]))
    A2 = np.vstack((np.sqrt(beta),alpha-xl2))
    J2 = np.vstack((A2[:,0:-1],A2[0,1:]))
    g1 = solve_banded((1,1),J1,en)
    g2 = solve_banded((1,1),J2,en)
    C = np.array(((1,-g1[-1]),(1,-g2[-1])))
    xl = np.array((xl1,xl2))
    ab = solve(C,xl)

    alphal = alpha
    alphal[-1] = ab[0]
    betal = beta
    betal[-1]=ab[1]
    x,w = __gauss__(alphal,betal)
    return x,w


def __rec_jacobi__(N,a,b):
    """ Generate the recursion coefficients alpha_k, beta_k

        P_{k+1}(x) = (x-alpha_k)*P_{k}(x) - beta_k P_{k-1}(x)

        for the Jacobi polynomials which are orthogonal on [-1,1]
        with respect to the weight w(x)=[(1-x)^a]*[(1+x)^b]

        Inputs:
        N - polynomial order
        a - weight parameter
        b - weight parameter

        Outputs:
        alpha - recursion coefficients
        beta - recursion coefficients

        Adapted from the MATLAB code by Dirk Laurie and Walter Gautschi
        http://www.cs.purdue.edu/archives/2002/wxg/codes/r_jacobi.m
    """

    from scipy.special import gamma

    nu = (b-a)/float(a+b+2)
    mu = 2**(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2)

    if N == 1:
        alpha = nu
        beta = mu
    else:
        n = np.arange(1.0,N)
        nab =  2*n+a+b
        alpha = np.hstack((nu,(b**2-a**2)/(nab*(nab+2))))
        n = n[1:]
        nab = nab[1:]
        B1 = 4*(a+1)*(b+1)/float((a+b+2)**2*(a+b+3))
        B = 4*(n+a)*(n+b)*n*(n+a+b)/(nab**2*(nab+1)*(nab-1))
        beta = np.hstack((mu,B1,B))

    return alpha, beta

def __polyval__(alpha,beta,x):
    """ Evaluate polynomials on x given the recursion coefficients
        alpha and beta """

    N = len(alpha)
    m = len(x)
    P = np.zeros((m,N+1))

    P[:,0] = 1
    P[:,1] = (x-alpha[0])*P[:,0]

    for k in range(1,N):
         P[:,k+1] = (x-alpha[k])*P[:,k] - beta[k]*P[:,k-1]

    return P

def __jacobi__(N,a,b,x):
    """ Compute the Jacobi polynomials which are orthogonal on [-1,1]
        with respect to the weight w(x)=[(1-x)^a]*[(1+x)^b]
        and evaluate them on the given grid up to P_N(x) """

    m = len(x)
    P = np.zeros((m,N+1))

    alpha, beta = __rec_jacobi__(N+1,a,b)

    P[:,0] = 1
    P[:,1] = (x-alpha[0])*P[:,0]

    for k in range(1,N):
         P[:,k+1] = (x-alpha[k])*P[:,k] - beta[k]*P[:,k-1]

    return P




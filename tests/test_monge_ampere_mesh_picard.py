# -*- coding: UTF-8 -*-
#! /usr/bin/python
from pigasus.utils.manager import context

# ...
try:
    from matplotlib import pyplot as plt
    PLOT=True
except ImportError:
    PLOT=False
# ...
from caid.cad_geometry import square
from caid.cad_geometry import circle
from caid.cad_geometry import quart_circle
from caid.cad_geometry import annulus
import numpy                as np
from time import time
import sys
import inspect
filename = inspect.getfile(inspect.currentframe()) # script filename (usually with path)
sys.stdout = open(filename.split('.py')[0]+'.txt', 'w')

# ... import picard from monge_ampere module
from pigasus.utils.load import load
monge_ampere    = load("monge_ampere")
picard          = monge_ampere.picard
# ...

abs = np.abs; sin = np.sin ; cos = np.cos ; exp = np.exp ; sqrt = np.sqrt
pi = np.pi; atan = np.arctan2 ; cosh = np.cosh
sech = lambda x: 1./cosh(x)

#-----------------------------------
try:
    nx = int(sys.argv[1])
except:
    nx = 15

try:
    ny = int(sys.argv[2])
except:
    ny = 15

try:
    px = int(sys.argv[3])
except:
    px = 2

try:
    py = int(sys.argv[4])
except:
    py = 2

geo   = square(n=[nx,ny], p=[px,py])
#geo   = circle(radius=1.,n=[nx,ny], p=[px,py])
#geo   = quart_circle(n=[nx,ny], p=[px,py])
#geo   = annulus(n=[nx,ny], p=[px,py])

#from caid.cad_geometry import cad_geometry as domain
#geo = domain("input/iter_inner.xml")
#-----------------------------------

#-----------------------------------
verbose = False
withTwoGrids = True

#    p_H = [ 5, 5]
#    p_h = [ 5, 5]

p_H = [ 3, 3]
p_h = [ 3, 3]

#    p_H = [ 2, 2]
#    p_h = [ 2, 2]

# TEST 1
#    # p = 5
#    rtol_H = 1.e-4
#    rtol_h = 1.e-8

# p = 3
#    rtol_H = 1.e-4
#    rtol_h = 1.e-6

# p = 2
#    rtol_H = 1.e-4
#    rtol_h = 1.e-4

# TEST 3
# p = 3, 5
rtol_H = 1.e-3
rtol_h = 1.e-3

#    # p = 2
##    rtol_H = 1.e-2
##    rtol_h = 1.e-2


rtol2_H = 1.e-6
#    rtol2_h = 1.e-6
rtol2_h = 1.e-9

maxiter_H = 40
maxiter_h = 40

n_H = [7,7]
nstage =  1
#-----------------------------------

#-----------------------------------
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

# ... test 2
#rho0 = lambda x,y : 1.
#C1   =  1.75484181939
#rho1 = lambda x,y : ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
# ...

# ... test 3
#rho0 = lambda x,y : 1.
#C1   = 0.285547502263
#rho1 = lambda x,y : (1. + 10*exp(-50*(y-0.5-0.25*sin(2*pi*x))**2))
# ...

# ...
#rho0 = lambda x,y : 1.
#t = 0.25
#C1 = 0.702563292151
#rho1 = lambda x,y : (1. + 5*exp(-50*abs((x-0.5-0.25*cos(2*pi*t))**2 \
#                    - (y-0.5-0.5 *sin(2*pi*t))**2 \
#                    - 0.01) ))
# ...

# ...
#rho0 = lambda x,y : 1.
#t = 1./3
#C1 = 0.831806957866
#rho1 = lambda x,y : ( 1. + 5*exp(-50*abs(y-0.5-0.25*sin(2*pi*x)*sin(2*pi*t))))
# ...

# ...
#rho0 = lambda x,y : 1.
#gamma = 5.
#lamb = 100.
#t = 0.75
#C1 = 0.832943327557
#x0 = t ; y0 = 0.2 + 0.5 * t ; x1 = 1. - t ; y1 = 0.8 - 0.5 * t
#u0 = lambda x,y : gamma * sech(lamb * ( x - x0 + y - y0 ))
#u1 = lambda x,y : gamma * sech(lamb * ( x - x1 + y - y1 ))
#rho1 = lambda x,y : ( 1. + u0(x,y) + u1(x,y))
# ...

# ... test7
#xc = 0.7 ; yc = 0.5
#C1 = 0.281648379406
#
#rho0 = lambda x,y : 1.
#r = lambda s,t : sqrt( (s-xc)**2 + (t-yc)**2 )
#theta = lambda s,t : atan(t-yc,s-xc)
#def rho1(s,t):
#    r_ = r(s,t) ;  t_ = theta(s,t)
#    val = C1 * (1. + 9./(1. + (10*r_*cos(t_-20*r_**2))**2) )
#    return val
# ...

## ...
# circle
#rho0 = lambda x,y : 1./pi
#C1   = 0.227475185932
#rho1 = lambda x,y : C1 * (1. + 5*exp(-25*abs((x-0.)**2+(y-0.)**2-0.4)))
## ...

# ...
# quart_circle
#rho0 = lambda x,y : 4./pi
#C1   = 2.91639889933
#rho1 = lambda x,y : C1 * ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
# ...

# ... annulus
#C0   = 0.424413181542
#rho0 = lambda x,y : C0 * 1.
#C1   = 0.733393862165
#rho1 = lambda x,y : C1 * ( 1. / (2. + cos(8*pi*sqrt((x-0.5)**2+(y-0.5)**2))))
# ...

# ...
#list_r  = np.genfromtxt('input/r.txt')
#list_ix = np.genfromtxt('input/ix.txt')
#list_iy = np.genfromtxt('input/iy.txt')
#
#list_kx = 2 * pi * np.asarray(list_ix)
#list_ky = 2 * pi * np.asarray(list_iy)
# ...

# ...
#C0   = 0.424413181542
#rho0 = lambda x,y : C0 * 1.
##C1   = 0.144432578196 # annulus, test 2
#C1   = 0.103412631611 # annulus, test 3
#def rho1(x,y):
#    window = exp(-10*abs((x-0.)**2+(y-0.)**2-0.8**2))
#    res = 0.
#    for (kx,ky,r) in zip(list_kx,list_ky,list_r):
##        res += r * (1+sin(kx*x)) * (1+sin(ky*y)) # test 2
#        res +=  (1+sin(x)) * (1+sin(ky*y)) # test 3
#    res *= window
#    res += 1.5
#    return res
# ...

# ...
#C0 = 0.0471623135665
#rho0 = lambda x,y : C0 * 1.
#C1 = 0.0223721108636 #ITER, test 2
##C1 = 0.0195621124256 #ITER, test 3
#def rho1(x,y):
#    window = exp(-10*abs(0.5*(x-6.)**2+0.15*(y-0.5)**2-1.0**2))
#    res = 0.
#    for (kx,ky,r) in zip(list_kx,list_ky,list_r):
#        res += r * (1+sin(kx*x)) * (1+sin(ky*y)) # test 2
##        res +=  (1+sin(x)) * (1+sin(ky*y)) # test 3
#    res *= window
#    res += 1.5
#    return res
# ...

# ...
n_h = []
for axis in range(0,2):
    n = n_H[axis]
    for i in range(0, nstage):
        n = 2*n+1
    n_h.append(n)

if withTwoGrids:
    print(">>>> coarse grid ", n_H, " with splines of degree ", p_H)
print(">>>> fine   grid ", n_h, " with splines of degree ", p_h)

if withTwoGrids:
    geo_H = square(n=n_H, p=p_H)
geo_h = square(n=n_h, p=p_h)
# ...

# ...
# values of gradu.n at the boundary
# ...
def func_g(x,y):
    return [x,y]
# ...

#-----------------------------------
# ...
# values of u at the boundary
# ...
bc_neumann={}
for data in geo_h.external_faces:
    patch_id = int(data[0]) ; face_id = int(data[1])
    bc_neumann[patch_id,face_id] = func_g
# ...
#-----------------------------------

# ...
tc = {}
tc['A'] = lambda x,y : [1., 0., 0., 1.]
tc['b'] = lambda x,y : [1.e-3]
tc['u'] = lambda x,y : [0.]
tc['f'] = lambda x,y : [0.]
tc['bc_neumann'] = bc_neumann
# ...

with context():

    #    PDE_H = picard(geometry=geo_H, testcase=tc)
    #    PDE_h = picard(geometry=geo_h, testcase=tc)
    if withTwoGrids:
        PDE_H = picard(geometry=geo_H, bc_neumann=bc_neumann)
    PDE_h = picard(geometry=geo_h, bc_neumann=bc_neumann)

    # ...
    print(">>> Solving using Picard <<<")
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
    c_rho = C0/C1
    # ...

    # ...
    if withTwoGrids:
        print("*****************************")
        tb = time()
        Errors_H, ErrorsH1_H = PDE_H.solve(  rho0, rho1, c_rho=None, u0=None \
                    , maxiter=maxiter_H, rtol=rtol_H, rtol2=rtol2_h, verbose=verbose)
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

    print("*****************************")
    tb = time()
    Errors_h, ErrorsH1_h = PDE_h.solve(  rho0, rho1, c_rho=None, u0=u0 \
                , maxiter=maxiter_h, rtol=rtol_h, rtol2=rtol2_h, verbose=verbose)
    te = time()
    print("Monge-Ampere eq. converges after ", len(Errors_h) \
            , " with final error ", Errors_h[-1] \
            , " with final H1-error ", ErrorsH1_h[-1])
    print("Elapsed time ", te-tb)
    print("*****************************")

    if withTwoGrids:
        uH = U_H.get()
    uh = U_h.get()

    if withTwoGrids:
        print("Error-coarse        ", np.abs(1.-PDE_H.norm(exact=PDE_H.Err_func)))
    print("Error-fine          ", np.abs(1.-PDE_h.norm(exact=PDE_h.Err_func)))

    if withTwoGrids:
        U_H.set(uH)
    U_h.set(uh)
    # ...

    # ...
    #    PDE_H.plotMesh(ntx=60, nty=60)
    # ...
    if PLOT:
        PDE_h.plotMesh(ntx=60, nty=60)
        plt.savefig(filename.split('.py')[0]+'.png', format='png')
        plt.clf()
    # ...

    np.savetxt("Errors.txt", np.asarray(Errors_h))

    if withTwoGrids:
        PDE_H.free()
    PDE_h.free()

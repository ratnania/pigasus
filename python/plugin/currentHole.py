# -*- coding: UTF-8 -*-
from pigasus.gallery.basicPDE import *
from pigasus.fem.field import *
from caid.cad_geometry import circle
import pylab as pl
import numpy as np

exp = np.exp ; log = np.log ; sqrt = np.sqrt

#-----------------------------------
nx              = 31
ny              = 31
px              = 2
py              = 2

geo = circle(n=[nx,ny], p=[px,py])
#-----------------------------------

#-----------------------------------
AllDirichlet    = True
nu              = 1.e-5
eta             = 1.e-5
jj1             = 0.2
jj2             = 0.266

dt              = 2.5e-2
niter           = 100
nfreq           = 10
t               = 0.
#-----------------------------------

#-----------------------------------
# Ew dictionary
tc_Ew = {}
tc_Ew['AllDirichlet'] = AllDirichlet
tc_Ew['b'] = lambda x,y : [1.]

# Iw dictionary
tc_Iw = {}
tc_Iw['AllDirichlet'] = AllDirichlet
tc_Iw['A'] = lambda x,y : [dt*nu, 0., 0., dt*nu]
tc_Iw['b'] = lambda x,y : [1.]

# Ephi dictionary
tc_Ephi = {}
tc_Ephi['AllDirichlet'] = AllDirichlet
tc_Ephi['b'] = lambda x,y : [-1.]

# Iphi dictionary
tc_Iphi = {}
tc_Iphi['AllDirichlet'] = AllDirichlet
tc_Iphi['A'] = lambda x,y : [1., 0., 0., 1.]

# Epsi dictionary
tc_Epsi = {}
tc_Epsi['AllDirichlet'] = AllDirichlet
tc_Epsi['b'] = lambda x,y : [1.]

# Ipsi dictionary
tc_Ipsi = {}
tc_Ipsi['AllDirichlet'] = AllDirichlet
tc_Ipsi['A'] = lambda x,y : [dt*eta, 0., 0., dt*eta]
tc_Ipsi['b'] = lambda x,y : [1.]

# Ej dictionary
tc_Ej = {}
tc_Ej['AllDirichlet'] = AllDirichlet
tc_Ej['A'] = lambda x,y : [-1., 0., 0., -1.]

# Ij dictionary
tc_Ij = {}
tc_Ij['AllDirichlet'] = AllDirichlet
tc_Ij['b'] = lambda x,y : [1.]
#-----------------------------------

# ------------------------------------------
# ...
Ew = basicPDE(geometry=geo, testcase=tc_Ew)
Iw = basicPDE(geometry=geo, testcase=tc_Iw)

wn   = Ew.rhs
wnew = Iw.unknown
# ...

# ...
Ephi = basicPDE(geometry=geo, testcase=tc_Ephi)
Iphi = basicPDE(geometry=geo, testcase=tc_Iphi)

phin   = Ephi.rhs
phinew = Iphi.unknown
# ...

# ...
Epsi = basicPDE(geometry=geo, testcase=tc_Epsi)
Ipsi = basicPDE(geometry=geo, testcase=tc_Ipsi)

psin   = Epsi.rhs
psinew = Ipsi.unknown
# ...

# ...
Ej = basicPDE(geometry=geo, testcase=tc_Ej)
Ij = basicPDE(geometry=geo, testcase=tc_Ij)

jn   = Ej.rhs
jnew = Ij.unknown
# ...
# ------------------------------------------

# ... set jc function
def func_jc ( x,y ) :
	R2 = x**2 + y**2
	lr_res = jj1 * ( 1.0 - R2**2 ) - jj2 * ( 1.0 - R2 )**8
	return [lr_res]

jn.set_func(func_jc)
# ...

# ... assembly all stationary operators

# ------------------------------------------
Ew.assembly()
Iw.assembly()

Ephi.assembly()
Iphi.assembly()

Epsi.assembly()
Ipsi.assembly()

Ej.assembly()
Ij.assembly()
# ------------------------------------------

# ... get jc as numpy array
jc = jn.get()
jc_values, = jn.evaluate(patch_id=0)
# ...

# ------------------------------------------
# redefine the right hand side function for Ew
def Fw(x,y):
    dxw  , dyw   = wn.grad(patch_id=0)
    dxj  , dyj   = jn.grad(patch_id=0)
    dxpsi, dypsi = psin.grad(patch_id=0)
    dxphi, dyphi = phin.grad(patch_id=0)

    v  = ( dxj * dypsi - dyj * dxpsi )
    v += ( dxphi * dyw - dyphi * dxw )
    v *= dt

    return [v]

Ew.rhs.set_func(Fw)
# ------------------------------------------

# ------------------------------------------
# redefine the right hand side function for Epsi
def Fpsi(x,y):
    dxpsi, dypsi = psin.grad(patch_id=0)
    dxphi, dyphi = phinew.grad(patch_id=0)

    v  = ( dxpsi * dyphi - dypsi * dxphi )
    v *= dt
    v -= eta * dt * jc_values

    return [v]

Epsi.rhs.set_func(Fpsi)
# ------------------------------------------

# ------------------------------------------
#           initialization
# ------------------------------------------
wn.reset()
wnew.reset()

phin.reset()
phinew.reset()

psin.reset()
psinew.reset()

jn.reset()
jnew.reset()

Mass        = Epsi.system
Stiffness   = Ej.system

j0          = Mass.solve(jc)
psi0        = Stiffness.solve(jc)

jn.set(j0)
psin.set(psi0)
# ------------------------------------------

# ------ Defining energy norms -------------
list_T  = [] # times array
list_Em = [] # magnetic energy
list_Ek = [] # kinetic energy
def magnetic_energy():
    lpr_Psi = psin.get()
    return np.dot ( lpr_Psi , Stiffness.dot(lpr_Psi) )

def kinetic_energy():
    lpr_Phi = phin.get()
    return np.dot ( lpr_Phi , Stiffness.dot(lpr_Phi) )
# ------------------------------------------

# ------------------------------------------
#           Fields Update
# ------------------------------------------
def update_j():
    jn.set(jnew)

def update_w():
    wn.set(wnew)

def update_phi():
    phin.set(phinew)

def update_psi():
    psin.set(psinew)

def update_fields():
    update_j()
    update_w()
    update_psi()
    update_phi()
# ------------------------------------------

# ------------------------------------------
#           LU factorization
# ------------------------------------------
Aw   = Iw.system.get().tocsc()
Apsi = Ipsi.system.get().tocsc()
Aphi = Iphi.system.get().tocsc()
Aj   = Ij.system.get().tocsc()

from scipy.sparse.linalg import splu

Aw_op   = splu(Aw)
Apsi_op = splu(Apsi)
Aphi_op = splu(Aphi)
Aj_op   = splu(Aj)
# ------------------------------------------

# ------------------------------------------
def solve(A_op, rhs, u):
    u.set(A_op.solve(rhs.get()))
# ------------------------------------------

# ------------------------------------------
#       Diagnostics routines
# ------------------------------------------
def diagnostics(i):
    if ( i % nfreq == 0 ):
        numdiag = i / nfreq
#        print("===		 diagnostic: "+str(numdiag)+"		===")

        lpr_Psi = psin.tomatrix(0)
        lpr_Phi = phin.tomatrix(0)
        lpr_W   = wn.tomatrix(0)
        lpr_J   = jn.tomatrix(0)

        np.savetxt("runs/J_"+str(numdiag)+".txt", lpr_J )
        np.savetxt("runs/Psi_"+str(numdiag)+".txt", lpr_Psi )
        np.savetxt("runs/Phi_"+str(numdiag)+".txt", lpr_Phi )
        np.savetxt("runs/w_"+str(numdiag)+".txt", lpr_W )

        # test si on ne diverge pas
        assert(lpr_J.max() > 1.e7)

def update_energies():
    list_T.append ( t )
    list_Em.append ( magnetic_energy ( ) )
    list_Ek.append ( kinetic_energy  ( ) )

def save_energies():
    np.savetxt( "runs/time.txt" , list_T )
    np.savetxt( "runs/Em.txt" , list_Em )
    np.savetxt( "runs/Ek.txt" , list_Ek )
# ------------------------------------------

update_energies()
# ------------------------------------------
for i in range(0, niter):
    t += dt

    # update w
    Ew.update()
    rhs = Ew.dot(wn) + Ew.rhs
    solve(Aw_op, rhs, wnew)     # Iw.solve(rhs)
    update_w()

    # update phi
    rhs = Ephi.dot(wnew)
    solve(Aphi_op, rhs, phinew) # Iphi.solve(rhs)
    update_phi()

    # update psi
    Epsi.update()
    rhs = Epsi.dot(psin) + Epsi.rhs
    solve(Apsi_op, rhs, psinew) # Ipsi.solve(rhs)
    update_psi()

    # update j
    rhs = Ej.dot(psinew)
    solve(Aj_op, rhs, jnew)     # Ij.solve(rhs)
    update_j()

    update_energies()
    diagnostics(i)

#Iphi.plot() ; pl.colorbar() ; pl.show()
save_energies()
# ------------------------------------------

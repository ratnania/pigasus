#########################################
nrb   = geo[0]

from caid.cad_geometry import cad_nurbs
C = np.zeros_like(nrb.points)
_C = np.genfromtxt("u.txt")
shape = list(nrb.shape)
C = np.zeros(shape+[3])
C[...,0] = _C
srf = cad_nurbs(nrb.knots, C, weights= nrb.weights)
#print srf.points

ntx = 80
nty = 80
#nty = 40
#ntx = 60
#nty = 60

tx = np.linspace(0., 1., ntx)
ty = np.linspace(0., 1., nty)

#tx = np.unique(nrb.knots[0])
#ty = np.unique(nrb.knots[1])

## ...
#P   = nrb.evaluate_deriv(tx,ty,nderiv=1)
#x   = P[0,:,:,0]
#xdu = P[1,:,:,0]
#xdv = P[2,:,:,0]
#
#y   = P[0,:,:,1]
#ydu = P[1,:,:,1]
#ydv = P[2,:,:,1]
#
#jac = xdu * ydv - xdv * ydu
## ...
#
## ...
#D   = srf.evaluate_deriv(tx,ty,nderiv=1)
#Udu = D[1,...,0]
#Udv = D[2,...,0]
#
#Udx =   ydv * Udu - ydu * Udv
#Udx /= jac
#Udy = - xdv * Udu + xdu * Udv
#Udy /= jac
## ...

# ...
P    = nrb.evaluate_deriv(tx,ty,nderiv=2)
x    = P[0,:,:,0]
xdu  = P[1,:,:,0]
xdv  = P[2,:,:,0]
xduu = P[3,:,:,0]
xduv = P[4,:,:,0]
xdvv = P[5,:,:,0]

y    = P[0,:,:,1]
ydu  = P[1,:,:,1]
ydv  = P[2,:,:,1]
yduu = P[3,:,:,1]
yduv = P[4,:,:,1]
ydvv = P[5,:,:,1]

jac = xdu * ydv - xdv * ydu
# ...

# ...
D    = srf.evaluate_deriv(tx,ty,nderiv=2)
Udu  = D[1,...,0]
Udv  = D[2,...,0]
Uduu = D[3,...,0]
Uduv = D[4,...,0]
Udvv = D[5,...,0]

Udx =   ydv * Udu - ydu * Udv
Udx /= jac
Udy = - xdv * Udu + xdu * Udv
Udy /= jac


C1 = Uduu - xduu * Udx - yduu * Udy
C2 = Uduv - xduv * Udx - yduv * Udy
C3 = Udvv - xdvv * Udx - ydvv * Udy
Udxx =   C1 * ydv**2    - 2 * C2 * ydu * ydv + C3 * ydu**2
Udxy = - C1 * xdv * ydv + C2 *(xdu * ydv + xdv * ydu) - C3 * xdu * ydu
Udyy =   C1 * xdv**2    - 2 * C2 * xdu * xdv + C3 * xdu**2
# ...


# ...
#P = srf(u=tx,v=ty)
#x = P[:,:,0]
#y = P[:,:,1]

fig = plt.figure()

Udx[:,0] = 0.
Udx[:,-1] = 1.
Udy[0,:] = 0.
Udy[-1,:] = 1.

for i,v in enumerate(ty):
#    phidx = Udu[:,i]
#    phidy = Udv[:,i]

    phidx = Udx[:,i]
    phidy = Udy[:,i]

    plt.plot(phidx, phidy, '-b')

for i,u in enumerate(tx):
#    phidx = Udu[i,:]
#    phidy = Udv[i,:]

    phidx = Udx[i,:]
    phidy = Udy[i,:]

    plt.plot(phidx, phidy, '-b')

plt.show()
#########################################

from caid.cad_geometry import square as domain
from caid.cad_geometry import cad_nurbs
import matplotlib.pyplot as plt
import numpy as np

pi = np.pi ; cos = np.cos ; sin = np.sin

n = 15
p = 2
geo = domain(n=[n,n], p=[p,p])
nrb = geo[0]

eps = 0.5 * (nrb.knots[0][p+1] - nrb.knots[0][0])
kx = 4 * pi
ky = 8 * pi

# ---------------------------------------------------
def genNewMesh(nrb):
    P = nrb.points

    x = P[:,:,0]
    y = P[:,:,1]

    shp = P.shape
    nx = shp[0] ; ny = shp[1]
    tx = np.linspace(0.,1.,nx)
    ty = np.linspace(0.,1.,ny)
    Tx,Ty = np.meshgrid(tx,ty)
    Tx = Tx.transpose()
    Ty = Ty.transpose()

    x[1:-1,1:-1] += eps*sin(kx*t*Tx[1:-1,1:-1])
    y[1:-1,1:-1] += eps*sin(ky*t*Ty[1:-1,1:-1])

    P[:,:,0] = x
    P[:,:,1] = y

    return cad_nurbs(nrb.knots, P, weights=nrb.weights)

def plotMesh(nrb):
    G = nrb.evalMesh(10)
    for d in G:
        plt.plot(d[:,0], d[:,1], '-k')
#    plt.show()
# ---------------------------------------------------

for t in np.linspace(0.,1.,10):
    nrb_new = genNewMesh(nrb)
    plotMesh(nrb_new)
    plt.savefig("mesh_t_"+str(t)+'.png', format='png')
    plt.clf()

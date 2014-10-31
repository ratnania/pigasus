from caid.cad_geometry import square as domain
from caid.cad_geometry import cad_geometry, cad_nurbs
import numpy as np
import sys

# shape of the input surface
nx_in = int(sys.argv[1])
ny_in = int(sys.argv[2])

# degrees of the input surface
px_in = int(sys.argv[3])
py_in = int(sys.argv[4])

# number of levels of the output surface
nstagex = int(sys.argv[5])
nstagey = int(sys.argv[6])

# degrees of the output surface
px = int(sys.argv[7])
py = int(sys.argv[8])

nx = nx_in ; ny = ny_in
for i in range(0, nstagex):
    nx = 2*nx+1
for i in range(0, nstagey):
    ny = 2*ny+1

geo = domain(n=[nx_in,ny_in], p=[px_in,py_in])
nrb   = geo[0]

C = np.zeros_like(nrb.points)
# import the solution
_C = np.genfromtxt("u.txt")
shape = list(nrb.shape)
C = np.zeros(shape+[3])
C[...,0] = _C
srf = cad_nurbs(nrb.knots, C, weights= nrb.weights)

geo_f = cad_geometry()
geo_f.append(srf)

geo_tmp = domain(n=[nx,ny], p=[px,py])
tx = [t for t in geo_tmp[0].knots[0] if t not in geo[0].knots[0]]
ty = [t for t in geo_tmp[0].knots[1] if t not in geo[0].knots[1]]

geo_f.refine(list_t = [tx,ty], list_p=[px-px_in, py-py_in])
u = geo_f[0].points[...,0]
import pylab as pl
pl.contourf(u) ; pl.colorbar() ; pl.show()
print(u.shape)
np.savetxt("u_ini.txt", u)

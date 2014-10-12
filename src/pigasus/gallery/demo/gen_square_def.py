import numpy as np
from caid.cad_geometry import square, cad_nurbs, cad_geometry
#from igakit.nurbs import NURBS

sqr = square(n=[0,0],p=[2,2])[0]
U,V = sqr.knots
C   = sqr.points

s = 1./np.sqrt(2)
weights         = np.ones((3,3))
weights[1,0]    = s
weights[0,1]    = s
weights[2,1]    = s
weights[1,2]    = s
srf = cad_nurbs([U,V], C, weights=weights)
geo = cad_geometry()
geo.append(srf)

try:
    nx, ny
except:
    nx = 3 ; ny = 5

tx = np.linspace(0.,1.,nx+2)[1:-1]
ty = np.linspace(0.,1.,ny+2)[1:-1]
geo.refine(0,list_t=[tx,ty])
geo._internal_faces = []
geo._external_faces = [[0,0],[0,1],[0,2],[0,3]]
geo._connectivity   = []

import sys
import numpy as np
import caid.cad_geometry as cg

nx = int(sys.argv[1])
ny = int(sys.argv[2])

px = int(sys.argv[3])
py = int(sys.argv[4])

mx = int(sys.argv[5])
my = int(sys.argv[6])

file_in = sys.argv[7]
file_out = sys.argv[8]

if (mx>px) or (my>py):
    raise("the multiplicity must be less than the final degree")

geo = cg.cad_geometry(file_in)

for i in range(0, geo.npatchs):
    list_p = geo[i].degree
    geo.refine(id=i,list_p=[px-list_p[0], py-list_p[1]])
    ux = np.linspace(0.,1.,nx+2)[1:-1]
    uy = np.linspace(0.,1.,ny+2)[1:-1]
    tx = []
    for x in ux:
        for j in range(0,mx):
            tx.append(x)
    ty = []
    for y in uy:
        for j in range(0,my):
            ty.append(y)

    tx = np.asarray(tx)
    ty = np.asarray(ty)

    geo.refine(id=i, list_t=[tx, ty])

geo.save(file_out)



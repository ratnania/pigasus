import sys
import numpy as np
import caid.cad_geometry as cg

nx = int(sys.argv[1])
px = int(sys.argv[2])

try:
    ls_file_in = sys.argv[3]
except:
    ls_file_in  = "domain_ini.xml"

try:
    ls_file_out = sys.argv[4]
except:
    ls_file_out  = "domain.xml"

geo = cg.cad_geometry(ls_file_in)

for i in range(0, geo.npatchs):
    list_p = geo[i].degree
    geo.refine(id=i,list_p=[px-list_p[0]])
    tx = np.linspace(0.,1.,nx+2)[1:-1]
    geo.refine(id=i, list_t=[tx])

geo.save(ls_file_out)



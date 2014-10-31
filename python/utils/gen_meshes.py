# -*- coding: UTF-8 -*-
import sys
import numpy as np
import caid.cad_geometry as cg

# ...
try:
    ls_file_in = sys.argv[1]
except:
    print("You must specify a domain description file.")
    sys.exit(0)
# ...

# ...
geo = cg.cad_geometry(ls_file_in)
dim = geo.dim
# ...

# ...
if dim==1:
    try:
        nin = [int(sys.argv[2])]
        pin = [int(sys.argv[3])]
    except:
        print("Incompatible data. The script will stop!")
        sys.exit(0)

    try:
        ls_file_out = sys.argv[4]
    except:
        ls_file_out = "domain.xml"

if dim==2:
    try:
        nin = [int(sys.argv[2]), int(sys.argv[3])]
        pin = [int(sys.argv[4]), int(sys.argv[5])]
    except:
        print("Incompatible data. The script will stop!")
        sys.exit(0)

    try:
        ls_file_out = sys.argv[6]
    except:
        ls_file_out = "domain.xml"

if dim==3:
    try:
        nin = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])]
        pin = [int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])]
    except:
        print("Incompatible data. The script will stop!")
        sys.exit(0)

    try:
        ls_file_out = sys.argv[8]
    except:
        ls_file_out = "domain.xml"
# ...

# ...
for i in range(0, geo.npatchs):
    nrb = geo[i]
    list_p = [p-pg for (p,pg) in zip(pin,nrb.degree)]
    # ...
    if min(list_p) < 0:
        print("specified degree is less than the original one")
        sys.exit(0)
    # ...
    geo.refine(id=i,list_p=list_p)
    list_t = []
    for (n,list_u) in zip(nin, nrb.knots):
        ub = list_u[0] ; ue = list_u[-1]
        t = np.linspace(ub,ue,n+2)[1:-1]
        list_t.append(t)
    geo.refine(id=i, list_t=list_t)
# ...

geo.save(ls_file_out)



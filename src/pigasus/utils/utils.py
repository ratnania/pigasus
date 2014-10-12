# -*- coding: UTF-8 -*-
import numpy as np

def hierarchical_geometries(geo, nlevel, domain):
    nh = geo[0].shape
    ph = geo[0].degree

    list_n = []
    for i in range(0,nlevel)[::-1]:
        nH = []
        for d in range(0, geo[0].dim):
            _n = (nh[d] - ph[d] - 1) / 2**i
            nH.append(_n)
        list_n.append(nH)

    list_geometry = []
    for n in list_n:
        geo = domain(n=n, p=ph)
        list_geometry.append(geo)

    return list_geometry

def process_diagnostics(sleep=None):
    def wait(sleep):
        import time
        print "Goes to sleep for "+str(sleep)+"s."
        time.sleep(sleep)  # Delay
    import os
    pid = os.getpid()
    print pid
    os.system("ps -p "+str(pid)+" -o %cpu,%mem,cmd")
    if sleep is not None:
        wait(sleep)

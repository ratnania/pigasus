# -*- coding: UTF-8 -*-
import numpy as np
import sys

# ------------------------------------------------------
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
# ------------------------------------------------------

# ------------------------------------------------------
def process_diagnostics(sleep=None):
    def wait(sleep):
        import time
        print("Goes to sleep for "+str(sleep)+"s.")
        time.sleep(sleep)  # Delay
    import os
    pid = os.getpid()
    print(pid)
    os.system("ps -p "+str(pid)+" -o %cpu,%mem,cmd")
    if sleep is not None:
        wait(sleep)
# ------------------------------------------------------

# ------------------------------------------------------
class progress():
    def __init__(self, n, title="", freq=2):
        self._n = n
        self._txt = " "
        self._title = " "+title
        self._freq = freq
#        sys.stdout.write(">> "+title+"\n")
        self.log(0)

    def finalize(self):
        self.log(self._n)
        self.reset()
        sys.stdout.write("\n")
#        sys.stdout.write("\ndone.\n")

    def log(self, i):
        if i % self._freq not in [0,1]:
            return
        self._i = i
        n = self._n
        self._txt += "."
#        txt = "\r[%d%%]"+self._txt
        txt = "\r[%d%%]"+self._title
        ratio = i * 100 / n
        sys.stdout.write(txt %ratio)
        sys.stdout.flush()

    def reset(self, n=None):
        self._txt = " "
        if n is not None:
            self._n = n
# ------------------------------------------------------

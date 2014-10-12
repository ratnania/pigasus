import matplotlib.pyplot    as plt
import numpy                as np

execfile('ex3.py')

xc = 0.65 ; yc = 0.
geo = geo_f.polarExtrude(xyzc=[xc,yc])
t = np.linspace(0.,1.,15+2)[1:-1]
geo.refine(list_t=[t,None])

geo.plotMesh()
plt.show()

import matplotlib.pyplot    as plt
import numpy                as np

execfile('ex5.py')

geo = geo_f.polarExtrude()
t = np.linspace(0.,1.,15+2)[1:-1]
geo.refine(list_t=[t,None])

geo.plotMesh()
plt.show()

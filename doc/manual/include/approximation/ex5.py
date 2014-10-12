import matplotlib.pyplot    as plt
import numpy                as np
from igakit.cad_geometry import cad_geometry
from igakit.cad_geometry import line as patch

from pigasus.fit.curfit import curfit
from pigasus.fit.curfit import compute_uk

cos = np.cos
sin = np.sin
sqrt = np.sqrt
pi = np.pi

#-----------------------------------
px          = 3
nx          = 63

alpha       = 1.

#method = "uniform"
method = "chord"
#method = "centripetal"
#-----------------------------------

#-----------------------------------
def MakeConstraint(cond, face=None, value=None):
    if cond.lower() == "closed":
        constraint = {}
        constraint['patch_id_m'] = 0
        constraint['face_m']     = 0
        constraint['patch_id_s'] = 0
        constraint['face_s']     = 1
    if cond.lower() == "c0":
        constraint = {}
        constraint['patch_id_m'] = 0
        constraint['face_m']     = face
        constraint['type']       = "C0"
        constraint['values']     = [value]
    if cond.lower() == "c1":
        constraint = {}
        constraint['patch_id_m'] = 0
        constraint['face_m']     = face
        constraint['type']       = "C1"
        constraint['values']     = [value]
    return constraint
#-----------------------------------

#-----------------------------------
# ...
tmin        = 0.
tmax        = 1.
npts        = 200
# ...

# ... import genPoints from grad_shafranov module
from pigasus.utils.load import load
grad_shafranov    = load("grad_shafranov")
genPoints         = grad_shafranov.genPoints
# ...

xyz = genPoints()

# ...
x,y = xyz
list_x = list(x) ; list_y = list(y)
# ...

# ...
list_Q = zip(list_x, list_y)
uk = compute_uk(list_Q, method=method)
U1       = []
U1      += list(uk)
list_xk  = []     ; list_yk  = []
list_xk += list_x ; list_yk += list_y

lists_uk = [U1]
lists_xk = [list_xk]
lists_yk = [list_yk]
# ...
#-----------------------------------

#-----------------------------------
# ....
constraints = []

# ... CLOSED
constraint = MakeConstraint("CLOSED")
constraints.append(constraint)
# ....
#-----------------------------------

#-----------------------------------
geo = patch(n=[nx], p=[px])
# ... make the patch periodic in 1 direction
geo._external_faces = [[0,0],[0,1]]

dict_con = {}
dict_con['original'] = [0,0]; dict_con['clone'] = [0,1]
geo._connectivity.append(dict_con)
#-----------------------------------

#-----------------------------------
fit = curfit(geometry=geo, constraints=constraints, alpha=alpha)
#-----------------------------------

#-----------------------------------
patch_id = 0
xk = lists_xk[patch_id]
yk = lists_yk[patch_id]

geo_f = fit.construct([xk, yk], uk=lists_uk)
#-----------------------------------

#-----------------------------------
srf = geo_f[0]

D = srf.evalMesh(10)

for d in D:
    plt.plot(d[:,0], d[:,1], '-k', label="constructed curve")

P = srf.points
nm = np.prod(np.asarray(srf.shape))
x = P[:,0].reshape(nm)
y = P[:,1].reshape(nm)
plt.plot(x,y,'Dr', label="Control points")

x = list_x
y = list_y
plt.plot(x, y, '.', color='#9999ff', label="Input data")

# Set x limits
plt.xlim(0.6, 1.6)
# Set y limits
plt.ylim(-0.6, 0.6)

plt.legend()
plt.show()
#-----------------------------------

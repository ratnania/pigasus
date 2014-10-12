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
nx          = 7

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
tmax        = 0.25
npts        = 50
# ...

# ... Quart Circle
f    = lambda s: [cos(2*pi*s), sin(2*pi*s)]
dsf  = lambda s: [-2*pi*sin(2*pi*s), 2*pi*cos(2*pi*s)]
d2sf = lambda s: [-(4*pi**2)*cos(2*pi*s), -(4*pi**2)*sin(2*pi*s)]
# ...

# ...
t = np.linspace(tmin,tmax, npts)
xyz = f(t)
# ...

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

# ... C0_COND
constraint = MakeConstraint("C0", 0, [1., 0.])
constraints.append(constraint)

constraint = MakeConstraint("C0", 1, [0., 1.])
constraints.append(constraint)

# ... C1_COND
constraint = MakeConstraint("C1", 0, [0., 1.])
constraints.append(constraint)

constraint = MakeConstraint("C1", 1, [1., 0.])
constraints.append(constraint)
# ....
#-----------------------------------

#-----------------------------------
geo = patch(n=[nx], p=[px])
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
plt.plot(x, y, 'o', color='#9999ff', label="Input data")

plt.legend()
plt.show()
#-----------------------------------

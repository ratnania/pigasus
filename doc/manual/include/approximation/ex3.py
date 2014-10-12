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
#method = "chord"
method = "centripetal"
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
npts        = 300
# ...

# ... X-point
_f1 = lambda s : [2*s,sqrt(1.-2*s)*(2*s) ]
_f2 = lambda s : [2.-2*s,-(2.-2*s)*sqrt(2*s-1.) ]
_df1 = lambda s : [2,2*sqrt(1.-2*s)-(2*s)/sqrt(1.-2*s) ]
_df2 = lambda s : [-2,2*sqrt(2*s-1.)-(2.-2*s)/sqrt(2*s-1) ]
def f(ls):
    x = [] ; y = []
    for s in ls:
        if s <= 0.5:
            [sx,sy] = _f1(s)
        else:
            [sx,sy] = _f2(s)
        x.append(sx) ; y.append(sy)
    return x, y

def dsf(ls):
    x = [] ; y = []
    for s in ls:
        if s <= 0.5:
            [sx,sy] = _df1(s)
        elif s > 0.5:
            [sx,sy] = _df2(s)
        x.append(sx) ; y.append(sy)
    return x, y
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
plt.plot(x, y, 'o', color='#9999ff', label="Input data")

# Set x limits
plt.xlim(-0.1, 1.1)
# Set y limits
plt.ylim(-0.5, 0.5)

plt.legend()
plt.show()
#-----------------------------------

# -*- coding: UTF-8 -*-
import numpy as np

n = 10

G0 = np.random.random((n,n))
G1 = np.random.random((n,n))

F0 = np.random.random((n,n))
F1 = np.random.random((n,n))

list_G = [G0,G1]
list_F = [F0,F1]

s = len(list_G)

# ...
G = np.zeros((s,s))
F = np.zeros((s,s))
for j in range(0,s):
    Gj = list_G[j]
    Fj = list_F[j]
    for i in range(0,s):
        Gi = list_G[i]
        Fi = list_F[i]
        G[i,j] = np.trace(Gi.transpose()*Gj)
        F[i,j] = np.trace(Fi.transpose()*Fj)
# ...

# ...
GF = G*F
trGF = np.trace(GF)
# ...

# ...
#def cost(x):
#    a = x[:s]
#    b = x[s:]
#    v = trGF
#    v += - 2 * np.dot(a,GF.dot(b))
#    v += np.dot(a,G.dot(a)) * np.dot(b,F.dot(b))
#    return v
# ...

# ...
def cost(x):
    a = x[:s]
    b = x[s:]
    # ...
    v0 = 0.
    for j in range(0,s):
        for i in range(0,s):
            v0 += G[i,j]*F[i,j]
    # ...

    # ...
    v1 = 0.
    for i in range(0,s):
        v1 += np.dot(G[i,:],a) * np.dot(F[i,:],b)
    # ...

    # ...
    v20 = 0.
    for j in range(0,s):
        for i in range(0,s):
            v20 += G[i,j] * a[i] * a[j]
    # ...

    # ...
    v21 = 0.
    for j in range(0,s):
        for i in range(0,s):
            v21 += F[i,j] * b[i] * b[j]
    # ...

    return v0 - 2 * v1 + v20 * v21
# ...

# ...
from scipy.linalg import kron
#from scipy.linalg import norm
from numpy.linalg import norm
def verification(a,b, list_G, list_F):
    G = list_G[0] ; F = list_F[0]
    R = kron(G,F)
    A = a[0]*G ; B = b[0]*F
    for i in range(1,s):
        G = list_G[i] ; F = list_F[i]
        R += kron(F,G)
        A += a[i]*G ; B += b[i]*F
    AoB = kron(A,B)
    print("norm Error ", norm(AoB-R, 'fro'))
# ...

# ...
from scipy.optimize import minimize
#x0 = np.zeros(2*s)
x0 = np.random.random(2*s)

#print "============"
#print "Initial ", cost(x0)
#res = minimize(cost, x0, method='CG', options={'gtol': 1e-6, 'disp': True})
#x = res.x
#a = x[:s]
#b = x[s:]
#print "a ", a
#print "b ", b
# ...

print("============")
print("Initial ", cost(x0))
a = x0[:s]
b = x0[s:]
verification(a,b, list_G, list_F)

method = 'CG'
#method = 'BFGS'
res = minimize(cost, x0, method=method, options={'gtol': 1e-12, 'disp': True})
x = res.x
a = x[:s]
b = x[s:]
# ...

# ...
verification(a,b, list_G, list_F)
# ...


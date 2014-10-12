import numpy as np
import pylab as pl

from numpy import pi, sin, cos

# ***************************************************
#
#                       TEST 102
#
# ***************************************************
eps = 1.e-1
m = 2 ; n = 1

f = lambda y :  sin( n*pi * ( 4*y*(1-y) ) )
df = lambda y : (-4*pi*n*y + 4*pi*n*(-y + 1))*cos(4*pi*n*y*(-y + 1))

u = lambda x,y : eps*sin(n*pi*y)*cos(m*pi*x) + f(y)

t = np.linspace(0.,1.,200)
x,y = np.meshgrid(t,t)

# ...
# plot of f
# ...
pl.plot(t, f(t))
pl.title("plot of $f(y) = \sin( n\pi ( 4 y (1-y) ) )$ ")
pl.show()
# ...

# ...
# plot of df
# ...
pl.plot(t, df(t))
pl.title("plot of $df(y) = -4 \pi n y + 4 \pi n (-y + 1)) \cos(4 \pi n y (-y + 1) $ ")
pl.show()
# ...

# ...
# plot of u
# ...
pl.contourf(x,y,u(x,y))
pl.colorbar()
pl.title("plot of $u (x,y) = \epsilon \sin(n \pi y) \cos(m \pi x) + \sin( n \pi ( 4 y (1-y) ) )$ ")
pl.show()
# ...

# ...
# plot of u - limit
# ...
eps = 0.
pl.contourf(x,y,u(x,y))
pl.colorbar()
pl.title("plot of $u (x,y) =  \sin( n \pi ( 4 y (1-y) ) )$ ")
pl.show()
# ...

# ***************************************************


# ***************************************************
#
#                       TEST 101
#
# ***************************************************
eps = 1.e-1
m = 2 ; n = 1

u = lambda x,y : sin(n*pi*y) + eps * cos(m*pi*x) * sin(n*pi*y)

t = np.linspace(0.,1.,201)
x,y = np.meshgrid(t,t)

# ...
# plot of u
# ...
pl.contourf(x,y,u(x,y))
pl.colorbar()
pl.title("plot of $u (x,y) = \sin(n \pi y) ( 1 + \epsilon \cos(m \pi x) )$ ")
pl.show()
# ...

# ...
# plot of u limit
# ...
eps = 0.
pl.contourf(x,y,u(x,y))
pl.colorbar()
pl.title("plot of $u (x,y) = \sin(n \pi y) $ ")
pl.show()
# ...

# ***************************************************


# ***************************************************
#
#                       TEST 110
#
# ***************************************************
eps = 1.e-4
alpha = 2.0
m = 2
n = 1

u = lambda x,y : eps*sin(n*pi*y)*cos(m*pi*x) + sin(alpha*(y**2 - y)*cos(m*pi*x) + n*pi*y)

t = np.linspace(0.,1.,201)
x,y = np.meshgrid(t,t)

# ...
# plot of u
# ...
pl.contourf(x,y,u(x,y))
pl.colorbar()
pl.title("plot of $u (x,y) = \epsilon \sin(n \pi y) \cos(m \pi x) + \sin( \\alpha (y^2 - y) \cos(m \pi x) + n \pi y)$ ")
pl.show()
# ...

# ...
# plot of u - limit
# ...
eps = 0.
pl.contourf(x,y,u(x,y))
pl.colorbar()
pl.title("plot of $u (x,y) = \sin( \\alpha (y^2 - y) \cos(m \pi x) + n \pi y)$ ")
pl.show()
# ...

# ***************************************************

import numpy as np

α = 1.0

def u_exact(x, y, t):
    return np.exp(-t)*np.exp(-x**2 - y**2)

def f_source(x, y, t):
    return (4*x**2 + 4*y**2 - 1)*np.exp(-t)*np.cos(2*x*y)

def initial_cond(x,y):
    return u_exact(x,y,0)

def bx0( y, t, x0):
    return u_exact(x0, y, t)

def bxf( y, t, xf):
    return u_exact(xf, y, t)

def by0( x, t, y0):
    return u_exact(x, y0, t)

def byf( x, t, yf):
    return u_exact(x, yf, t)

x0 = -4.0
xf = 4.0
y0 = -4.0
yf = 4.0

t0 = 0.0
tf = 5000
Nx = 40
Ny = 41
NtimeSteps = 50


#
dt = tf/NtimeSteps
dx = (xf - x0)/Nx
dy = (yf - y0)/Ny

x = np.linspace(x0, xf, Nx+1)
y = np.linspace(y0, yf, Ny+1)
t = np.linspace(t0, tf, NtimeSteps+1)

X, Y = np.meshgrid(x, y)

u = np.zeros((Nx+1,Ny+1))

# Initial condition
u[:,:] = initial_cond(X, Y).T


import matplotlib.pyplot as plt

fig = plt.figure()

ax1 = fig.add_subplot(121, projection="3d")
t = 0.0
ax1.plot_surface(X, Y, u_exact(X, Y, t), cmap="coolwarm" )
str_title = "t = {:10.5f}".format(t)
ax1.set_title(str_title)
ax1.set_zlim(0.0, 0.9)
ax1.set_xlabel("x")
ax1.set_ylabel("y")


ax2 = fig.add_subplot(122, projection="3d")
t = 0.5
ax2.plot_surface(X, Y, u_exact(X, Y, t), cmap="coolwarm" )
str_title = "t = {:10.5f}".format(t)
ax2.set_title(str_title)
ax2.set_zlim(0.0, 0.9)
ax2.set_xlabel("x")
ax2.set_ylabel("y")

plt.show()


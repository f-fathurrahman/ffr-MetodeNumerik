import numpy as np
import matplotlib.pyplot as plt

L = 1.5
T = 1.0
alpha = 0.5

def u_exact(t, x):
    return 5*x*t*(L-x)

def initial_cond(x):
    return u_exact(0, x)

def source_term(t, x):
    return 10*alpha*t + 5*x*(L - x)

Nx = 50
x = np.linspace(0.0, L, Nx+1)
dx = x[1] - x[0]
print("dx = ", dx)

Nt = 100
t = np.linspace(0.0, T, Nt+1)
dt = t[1] - t[0]
print("dt = ", dt)

F = alpha * dt / dx**2
print("F = ", F)
if F > 0.5:
    print("WARNING: solution is not stable")


#plt.clf()
#plt.plot(x, initial_cond(x), label="initial")
#plt.savefig("IMG_diff_0.png", dpi=150)

u = np.zeros((Nt+1, Nx+1))

u[0,:] = initial_cond(x)
for n in range(0,Nt):
    u[n,0] = 0.0 # syarat batas pada x=0
    u[n,Nx] = 0.0 # syarat batas pada x=L
    for i in range(1,Nx):
        fni = source_term(t[n], x[i])
        u[n+1,i] = u[n,i] + F*( u[n,i+1] - 2*u[n,i] + u[n,i-1] ) + dt *  fni

    plt.clf()
    plt.plot(x, u[n,:])
    plt.title("t = " + str(t[n]))
    plt.ylim(0.0, 0.3)
    plt.savefig("IMG_diff1d_explicit" + str(n) + ".png", dpi=150)

    print("n = " + str(n) + " is done")

plt.clf()
plt.plot(x, u[Nt,:], label="numerical")
plt.plot(x, u_exact(t[Nt],x), label="exact sol")
plt.title("t = " + str(t[Nt]))
plt.ylim(0.0, 0.3)
plt.legend()
plt.savefig("IMG_diff1d_explicit_COMPARE_" + str(n) + ".png", dpi=150)

du = u_exact(t[Nt],x) - u[Nt,:]
norm_du = np.linalg.norm(du)
print("norm_du = ", norm_du)
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

Nt = 50
t = np.linspace(0.0, T, Nt+1)
dt = t[1] - t[0]
print("dt = ", dt)

F = alpha * dt / dx**2
print("F = ", F)

u = np.zeros((Nt+1, Nx+1))
# Set syarat awal
u[0,:] = 0.0

# Build matrix A
A = np.zeros( (Nx-1,Nx-1) )

A[0,0] = 1 + 2*F
A[0,1] = -F
for i in range(1,Nx-2):
    A[i,i] = 1 + 2*F
    A[i,i+1] = -F
    A[i,i-1] = -F
A[Nx-2,Nx-2] = 1 + 2*F
A[Nx-2,Nx-3] = -F

# time loop
b = np.zeros(Nx-1)
for n in range(1,Nt+1):
    # Vektor RHS
    b[:] = 0.0
    for i in range(Nx-1):
        b[i] = u[n-1,i+1] + dt*source_term(t[n], x[i+1])
    # Selesaikan sistem persamaan linear
    u[n,1:Nx] = np.linalg.solve(A, b)

    plt.clf()
    plt.plot(x, u[n,:])
    plt.title("t = " + str(t[n]))
    plt.ylim(0.0, 5.0)
    plt.savefig("IMG_diff1d_implicit_" + str(n) + ".png", dpi=150)

    print("n = " + str(n) + " is done")

plt.clf()
plt.plot(x, u[Nt,:], label="numerical")
plt.plot(x, u_exact(t[Nt],x), label="exact sol")
plt.title("t = " + str(t[Nt]))
plt.ylim(0.0, 5.0)
plt.legend()
plt.savefig("IMG_diff1d_implicit_COMPARE_" + str(n) + ".png", dpi=150)

du = u_exact(t[Nt],x) - u[Nt,:]
norm_du = np.linalg.norm(du)
print("norm_du = ", norm_du)


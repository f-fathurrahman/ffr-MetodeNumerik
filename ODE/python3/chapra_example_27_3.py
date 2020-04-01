# Finite-difference method for linear BVP

# d2T/dx2 + h'(T_a - T) = 0
# Boundary condition:
#   T(0) = 40
#   T(10) = 200

import numpy as np

h = 0.01
T_a = 20.0

x0 = 0.0
T0 = 40.0 # Boundary condition, T(0) = 40
xf = 10.0
Tf = 200.0 # Boundary condition, T(10) = 200
Δx = 2.0 # segment length (or step size in shooting method)
Nstep = int( (xf-x0)/Δx )
Npoints = Nstep + 1

T = np.zeros(Npoints)
T[0] = T0
T[-1] = Tf

# Finite-difference operator of second derivative matrix
# In general we should use sparse matrix. However, because
# The size is rather small, we use full (dense) matrix
# Please refer to the left-hand-side of Eq. 27.3 for the matrix elements.
Npointsm2 = Npoints-2 # Number of interior points
d2dx2 =  np.zeros((Npointsm2,Npointsm2))
for i in range(Npointsm2):
    d2dx2[i,i] = 2 + h*Δx**2
    if i != 0:
        d2dx2[i-1,i] = -1.0
    if i != (Npointsm2-1):
        d2dx2[i+1,i] = -1.0
# Display the matrix
print("FD representation of second-derivative operator:")
print(d2dx2)

# The vector represented by the right hand side of Eq. 27.3
f = np.zeros(Npointsm2)
for i in range(1,Npointsm2-1):
    f[i] = h*Δx**2*T_a
# From the left BC
f[0] = h*Δx**2*T_a + T0
# From the right BC 
f[-1] = h*Δx**2*T_a + Tf
# Display
print(f)

# Solve the linear equations
T[1:Npoints-1] = np.linalg.solve(d2dx2,f)
x = np.zeros(Npoints)
for i in range(Npoints):
    x[i] = x0 + i*Δx
    print("%f %f" % (x[i], T[i]))

import matplotlib.pyplot as plt
plt.clf()
plt.plot(x, T, marker="o", label="Temperature")
plt.xlabel("x")
plt.ylabel("T")
plt.legend()
plt.savefig("IMG_example_27_3.png", dpi=150)
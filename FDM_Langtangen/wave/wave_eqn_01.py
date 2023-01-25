import numpy as np
import matplotlib.pyplot as plt

# Definisi domain
T = 10.0 # waktu
L = 10.0 # panjang

# Parameter diskritisasi
Nx = 50
Nt = 100

x = np.linspace(0.0, L, Nx+1) # variabel independen
t = np.linspace(0.0, T, Nt+1) # variabel independen

dx = x[1] - x[0]
dt = t[1] - t[0]

c = 1.5
print("dx = ", dx)
print("dt = ", dt)

# Bilangan Courant
C = c * dt / dx
print("C = ", C)
if C > 1:
    print("WARNING: iterasi tidak stabil")

def initial_func(x):
    return np.sin(4*np.pi*x/L)

def segitiga(x):
    if x <= 5.0:
        return x/5
    else:
        return 2.0 - x/5

# Alokasi array untuk solusi
u = np.zeros( (Nt+1,Nx+1) )

# Aplikasi syarat awal
u[0,:] = initial_func(x)

#for i in range(Nx):
#    u[0,i] = segitiga(x[i])

# Plot untuk t=0
plt.clf()
plt.plot(x, u[0,:])
plt.ylim(-1.1, 1.1)
plt.title("t = " + str(t[0]))
plt.savefig("IMG_wave_t_" + str(0) + ".png", dpi=150)

#exit()

# Step pertama
# Aplikasi syarat batas (kiri dan kanan)
u[1,0] = 0.0 # kiri
u[1,Nx] = 0.0 # kanan
for i in range(1,Nx):
    u[1,i] = u[0,i] + 0.5*C**2 * ( u[0,i+1] - 2*u[0,i] + u[0,i-1] )


# Plot untuk t=t[n]
n = 1 # indeks waktu
plt.clf()
plt.plot(x, u[n,:])
plt.ylim(-1.1, 1.1)
plt.title("t = " + str(t[n]))
plt.savefig("IMG_wave_t_" + str(n) + ".png", dpi=150)


for n in range(1,Nt):
    u[n+1,0] = 0.0 # kiri
    u[n+1,Nx] = 0.0 # kanan
    for i in range(1,Nx):
        u[n+1,i] = 2*u[n,i] - u[n-1,i] + C**2 * ( u[n,i+1] - 2*u[n,i] + u[n,i-1] )

    plt.clf()
    plt.plot(x, u[n+1,:])
    plt.ylim(-1.1, 1.1)
    plt.title("t = " + str(t[n+1]))
    plt.savefig("IMG_wave_t_" + str(n+1) + ".png", dpi=150)

    print("n = ", n, " sudah selesai")


print("Program selesai")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

Xmax = 401
Ymax = 101
Zmaz = 100

ε = 4
dd = 0.5
Ex = np.zeros(Xmax)
Hy = np.zeros(Xmax)
β = np.zeros(Xmax)

for i in range(401):
    if i < 201:
        β[i] = dd
    else:
        β[i] = dd/ε

z = np.arange(201)
zx = np.arange(1,Xmax-1)
xs = np.arange(1,Xmax-1)
Ex[:201] = 0.5*np.sin(2*np.pi*z/100.0)
Hy[:201] = 0.5*np.sin(2*np.pi*z/100.0)

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(1,Xmax-1), ylim=(-1.5,1.5))
ax.grid()
line, = ax.plot(xs, Ex[1:Xmax-1], lw=2)

def animate(dum):
    for x in range(1,Xmax-1):
        Ex[x] = Ex[x] + β[x]*(Hy[x-1] - Hy[x])
        Hy[x] = Hy[x] + dd*(Ex[x] - Ex[x+1])
    line.set_data(xs, Ex[1:Xmax-1])
    return line, 

plt.title("Title")
plt.xlabel("z")
plt.ylabel("Ex")
p = plt.axvline(x=200, color="r") # Vertical line separator
ani = animation.FuncAnimation(fig , animate,1 , blit=True)
plt.show()
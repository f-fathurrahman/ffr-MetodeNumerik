import matplotlib.pyplot as plt
import numpy as np

def my_func(c):
    g = 9.81
    m = 68.1
    v = 40.0
    t = 10.0
    return g*m/c*(1 - np.exp(-(c/m)*t)) - v

c = np.linspace(13.0, 16.0, 200)
f = my_func(c)

xl = 14.0
xu = 15.0
xr = (xl + xu)/2

fl = my_func(xl)
fu = my_func(xu)
fr = my_func(xr)

print("Iterasi 1: fr = ", fr)

plt.clf()
plt.plot(c, f)
plt.grid(True)
plt.plot([xl], [fl], marker="o", label="xl", color="red")
plt.plot([xu], [fu], marker="o", label="xu", color="green")
plt.plot([xr], [fr], marker="o", label="xr", color="magenta")
plt.legend()
plt.title("Iterasi 1")
plt.savefig("IMG_01.png", dpi=150)
plt.show()


# Iterasi 2

xl = xr
xu = 15.0
xr = (xl + xu)/2

fl = my_func(xl)
fu = my_func(xu)
fr = my_func(xr)

print("Iterasi 2: fr = ", fr)

plt.clf()
plt.plot(c, f)
plt.grid(True)
plt.plot([xl], [fl], marker="o", label="xl", color="red")
plt.plot([xu], [fu], marker="o", label="xu", color="green")
plt.plot([xr], [fr], marker="o", label="xr", color="magenta")
plt.legend()
plt.savefig("IMG_02.png", dpi=150)
plt.title("Iterasi 2")
plt.show()


# Iterasi 3

xl = xr
xu = 15.0
xr = (xl + xu)/2

fl = my_func(xl)
fu = my_func(xu)
fr = my_func(xr)

print("Iterasi 3: fr = ", fr)

plt.clf()
plt.plot(c, f)
plt.grid(True)
plt.plot([xl], [fl], marker="o", label="xl", color="red")
plt.plot([xu], [fu], marker="o", label="xu", color="green")
plt.plot([xr], [fr], marker="o", label="xr", color="magenta")
plt.legend()
plt.title("Iterasi 3")
plt.savefig("IMG_03.png", dpi=150)
plt.show()



# Iterasi 4

xu = xr
xr = (xl + xu)/2

fl = my_func(xl)
fu = my_func(xu)
fr = my_func(xr)

print("Iterasi 4: fr = ", fr)

plt.clf()
plt.plot(c, f)
plt.grid(True)
plt.plot([xl], [fl], marker="o", label="xl", color="red")
plt.plot([xu], [fu], marker="o", label="xu", color="green")
plt.plot([xr], [fr], marker="o", label="xr", color="magenta")
plt.legend()
plt.title("Iterasi 4")
plt.savefig("IMG_04.png", dpi=150)
plt.show()


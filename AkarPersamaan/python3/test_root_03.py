import math
from search_root_interval import *
from bisection import *
from regula_falsi import *
from ridder import *
import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return x - np.tan(x)

a = -1.0
b = 20.0
dx = 0.01

plt.clf()
xplot = np.linspace(a,b,2000)
plt.plot(xplot, f(xplot))
plt.ylim([-1.0,1.0])
plt.savefig("TEMP_plot_03.png", dpi=300)

itrial = 0
while True:
    
    x1, x2 = search_root_interval(f,a,b,dx,verbose=True)
    
    if x1 != None:

        itrial = itrial + 1

        plt.clf()
        xplot = np.linspace(x1,x2,500)
        plt.plot(xplot, f(xplot))
        plt.savefig("TEMP_plot_" + str(itrial) + ".png", dpi=300)

        a = x2
        print("Check f(x1)f(x2) = ", f(x1)*f(x2))
        #root, _ = bisection(f, x1, x2, verbose=False)
        root, _ = ridder(f, x1, x2)
        if root != None:
            print("Root is found: %18.10f" % root)
            print("f(x) = ", f(root))
    else:
        print("\nDone")
        break


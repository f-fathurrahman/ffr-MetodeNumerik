from numpy import *
import matplotlib.pyplot as plt
import time
plt.ion()
x = linspace(0, 2*pi, 100)
line, = plt.plot(x, sin(x))
for b in linspace(0 , 4*pi , 40) :
    line.set_ydata(sin(x + b))
    plt.draw()
    time.sleep(0.1)

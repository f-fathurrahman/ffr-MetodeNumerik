from vpython import *
import numpy as np

Xm = 201
Ym = 100
Zm = 100
ts = 2
β = 0.01
Ex = np.zeros((Xm,ts))
Hy = np.zeros((Xm,ts))

scene = canvas(x=0, y=0, width=800, height=500,
    title="E: cyan, H: red, PBC", forward=vector(-0.6,-0.5,-1))

Eplot = curve(x=list(range(0,Xm)), color=color.cyan, radius=1.5, display=scene)

Hplot = curve(x=list(range(0,Xm)), color=color.red, radius=1.5, display=scene)

vplane = curve(vector(-Xm,Ym,0), vector(Xm,Ym,0),
    vector(Xm,-Ym,0),vector(-Xm,-Ym,0), vector(-Xm,Ym,0),
    color=color.cyan)

zaxis = curve(vector(-Xm,0,0), vector(Xm,0,0), color=color.magenta)

hplane = curve( vector(-Xm,0,Zm), vector(Xm,0,Zm),
    vector(Xm,0,-Zm), vector(-Xm,0,-Zm), vector(-Xm,0,Zm),
    color=color.magenta )

ball1 = sphere(pos=vector(Xm+30,0,0), color=color.black, radius=2)

ExLabel1 = label(text="Ex", pos=vector(-Xm-10,50,0), box=0)
ExLabel1 = label(text="Ex", pos=vector(Xm+10,50,0), box=0)
HyLabel = label(text="Hy", pos=vector(-Xm-10,0,50), box=0)
zLabel = label(text="Z", pos=vector(Xm+10,0,0), box=0)

def plot_fields():
    z = arange(Xm)
    for i in range(Xm):
        Eplot.x[i] = 2*z[i] - Xm
        Eplot.y[i] = 800*Ex[z[i],0]
        Hplot.x[i] = 2*z[i] - Xm
        Hplot.z[i] = 800*Hy[z[i],0]

z = arange(Xm)
for i in range(Xm):
    Ex[i,0] = 0.1*sin(2*pi*z[i]/100.0)
    Hy[i,0] = 0.1*sin(2*pi*z[i]/100.0)

print(Eplot.x)
plot_fields()
print(Eplot.x)


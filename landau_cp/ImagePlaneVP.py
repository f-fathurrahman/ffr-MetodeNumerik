import vpython as vpy

scene = vpy.canvas(width=500, height=500, range=100,
    title="E of Charge Left of Plane")
plane = vpy.box(pos=vpy.vector(0,0,0), length=2, height=130, width=130,
    color=vpy.vector(0.9,0.9,0.9), opacity=0.5)
gridpts = vpy.points(radius=0.01, color=vpy.color.cyan)
plusChg = vpy.sphere(radius=5, color=vpy.color.red, pos=vpy.vector(40,0,0))
negChg = vpy.sphere(radius=5, color=vpy.color.red, pos=vpy.vector(-40,0,0))

def grid3d():
    for z in range(-60,80,20):
        for y in range(-60,80,20):
            for x in range(-60,80,20):
                gridpts.append(pos=vpy.vector(x,y,z))

def electricF():
    for y in range(-60,80,20):
        for z in range(-60,80,20):
            for x in range(-60,80,20):
                r = vpy.vector(x,y,z)
                xx =vpy.vector(40.0,0,0)
                d = vpy.vector(r - xx)
                dm = vpy.vector(r + xx)
                dd = vpy.mag(d)
                ddp = vpy.mag(dm)
                if(x == 40 and y==0 and z==0):
                    continue
                if ddp != 0:
                    E1 = d/dd**3
                    E2 = -dm/ddp**3
                    E = E1 + E2
                    elecF = vpy.arrow(pos=r, color=vpy.color.orange)
                    elecF.axis = 10*E/vpy.mag(E)

grid3d()
electricF()


import numpy as np
import pyvista as pv

mesh = pv.read("dumbbell.vtk")

points = mesh.points

# Sphere parameters (same as Gmsh)
r = 1.0
L = 3.0
x1 = -L/2
x2 =  L/2

# Distance from sphere centers
d1 = np.linalg.norm(points - np.array([x1,0,0]), axis=1)
d2 = np.linalg.norm(points - np.array([x2,0,0]), axis=1)

# Points strictly inside spheres (not boundary)
inside = (d1 < r - 1e-3) | (d2 < r - 1e-3)

inside_points = points[inside]

plotter = pv.Plotter()

plotter.add_mesh(mesh, opacity=0.2)
plotter.add_points(inside_points, color="red", point_size=8)

plotter.show()

print("Number of internal points:", inside_points.shape[0])
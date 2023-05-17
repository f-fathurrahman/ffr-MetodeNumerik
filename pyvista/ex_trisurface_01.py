# https://docs.pyvista.org/version/stable/examples/00-load/create-tri-surface.html

import numpy as np
import pyvista as pv

# Define a simple Gaussian surface
n = 20
x = np.linspace(-200, 200, num=n) + np.random.uniform(-5, 5, size=n)
y = np.linspace(-200, 200, num=n) + np.random.uniform(-5, 5, size=n)
xx, yy = np.meshgrid(x, y)
A, b = 100, 100
zz = A * np.exp(-0.5 * ((xx / b) ** 2.0 + (yy / b) ** 2.0))

# Get the points as a 2D NumPy array (N by 3)
points = np.c_[xx.reshape(-1), yy.reshape(-1), zz.reshape(-1)]
print(points[0:5,:])

# simply pass the numpy points to the PolyData constructor
cloud = pv.PolyData(points)
cloud.plot(point_size=15)

surf = cloud.delaunay_2d()
surf.plot(show_edges=True)

#
# Masked Triangulations
#
x = np.arange(10, dtype=float)
xx, yy, zz = np.meshgrid(x, x, [0])
points = np.column_stack((xx.ravel(order="F"), yy.ravel(order="F"), zz.ravel(order="F")))

# Perturb the points
points[:, 0] += np.random.rand(len(points)) * 0.3
points[:, 1] += np.random.rand(len(points)) * 0.3

# Create the point cloud mesh to triangulate from the coordinates
cloud = pv.PolyData(points)
cloud


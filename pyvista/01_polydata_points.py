import pyvista as pv
import numpy as np

#points = np.random.random((1000, 3))
#pc = pv.PolyData(points)
#pc.plot(scalars=points[:, 2], point_size=5.0, cmap='jet')

points = np.zeros((4,3))
points[0,:] = [-1.0, 1.0, 0.0]
points[1,:] = [-3.0, 1.0, 0.0]
points[2,:] = [-1.0, 2.0, 0.0]
points[3,:] = [-4.0, 1.0, 0.0]
pc = pv.PolyData(points)
pc.plot(color="tan", point_size=10)
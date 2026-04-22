import pyvista as pv

mesh = pv.read("dumbbell.vtk")

plotter = pv.Plotter()

# Show mesh (semi-transparent helps)
plotter.add_mesh(mesh, opacity=0.3, show_edges=True)

# Show ALL points
plotter.add_points(mesh.points, render_points_as_spheres=True, point_size=6)

plotter.show()
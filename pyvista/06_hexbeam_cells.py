import pyvista as pv
from pyvista import examples

mesh = examples.load_hexbeam()

pl = pv.Plotter()
pl.add_mesh(mesh, show_edges=True, color="white", opacity=0.5)
pl.add_points(mesh.points, color="red", point_size=20)

#single_cell = mesh.extract_cells(mesh.n_cells - 1)
single_cell = mesh.extract_cells(0)
pl.add_mesh(single_cell, color="green", edge_color="blue",
            line_width=5, show_edges=True)

pl.camera_position = [(6.20, 3.00, 7.50),
                      (0.16, 0.13, 2.65),
                      (-0.28, 0.94, -0.21)]
pl.show()

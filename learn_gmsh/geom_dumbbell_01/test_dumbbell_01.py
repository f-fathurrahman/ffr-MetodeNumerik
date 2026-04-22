import pyvista as pv

# Load mesh (vtk or msh converted to vtk)
mesh = pv.read("dumbbell.vtk")

# --------------------------------------------------
# STEP 1: Extract ONLY outer surface
# --------------------------------------------------
surface = mesh.extract_surface(algorithm="dataset_surface")

# --------------------------------------------------
# STEP 2: Clean (remove duplicate/internal faces)
# --------------------------------------------------
surface = surface.clean()

# --------------------------------------------------
# STEP 3: Optional smoothing (nice visualization)
# --------------------------------------------------
surface = surface.smooth(n_iter=30)

# --------------------------------------------------
# Plot
# --------------------------------------------------
plotter = pv.Plotter()
plotter.add_mesh(surface, show_edges=True)
plotter.show()


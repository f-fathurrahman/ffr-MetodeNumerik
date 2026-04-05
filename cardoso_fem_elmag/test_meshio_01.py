import meshio

mesh = meshio.read("square2.msh")
print("Nodes:", mesh.points)
print("Triangles:", mesh.cells_dict["triangle"])
print("1D Groups:", mesh.cell_data_dict["gmsh:physical"]["line"])
print("2D Groups:", mesh.cell_data_dict["gmsh:physical"]["triangle"])
print("Group Names:")
for name, (tag, dim) in mesh.field_data.items():
    print(f"{name}: ID {tag}, dim {dim}")
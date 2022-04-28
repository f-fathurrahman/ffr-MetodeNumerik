import MAT

meshfile = "../TEMP_octave/mesh_cylinder.mat"
vars = MAT.matread(meshfile)
g_coord = vars["g_coord"]
g_num = convert(Matrix{Int64}, vars["g_num"])

using WriteVTK: MeshCell, VTKCellTypes, vtk_grid
points = g_coord
Ncells = size(g_num,2)
cells = Vector{MeshCell}(undef,Ncells)
for i in 1:Ncells
    @views cells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, g_num[:,i])
end



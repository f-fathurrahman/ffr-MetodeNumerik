using WriteVTK

x = 0.0:0.1:1.0
y = 0.0:0.2:1.0
z = -1:0.05:1.0

vtk_grid("TEMP_01", x, y, z) do vtk
    vtk["temperature"] = rand(length(x), length(y), length(z))
end
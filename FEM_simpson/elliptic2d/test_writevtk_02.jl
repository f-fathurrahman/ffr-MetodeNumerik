# Rectilinear data

using WriteVTK

x = [0.0, 0.1, 0.5, 1.3]
y = sort(rand(8))
z = [-cospi(i/10) for i = 0:10]

vtk_grid("TEMP_02", x, y, z) do vtk
    vtk["temperature"] = rand(length(x), length(y), length(z))
end
using WriteVTK

x = -5.0:0.1:5.0
y = -10.0:0.1:10.0
z = -5.0:0.1:5.0

Nx = length(x)
Ny = length(y)
Nz = length(z)

T = zeros(Nx,Ny,Nz)
for k in 1:Nx, j in 1:Ny, i in 1:Nx
    #T[i,j,k] = sin(x[i] * y[j]) * sin(y[j] * z[k]) * cos(z[k] * x[i])
    r2 = x[i]^2 + y[j]^2 + z[k]^2
    T[i,j,k] = exp(-0.1*r2)
end

vtk_grid("TEMP_03", x, y, z) do vtk
    vtk["temperature"] = T
end


x = linspace( -8.0, 8.0, 41 )';
y = linspace( -8.0, 8.0, 41 )';
[X, Y] = meshgrid(x,y)
Z = sqrt(X.^2 + Y.^2)

plot3d(x, y, Z)

if getscilabmode() ~= "STD"
  quit()
end

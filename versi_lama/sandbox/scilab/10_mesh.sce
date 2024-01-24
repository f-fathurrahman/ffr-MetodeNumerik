[X,Y] = meshgrid(-5:0.5:5, -5:0.5:5)
Z = sin(X) + cos(Y)*cos(X)

clf()
mesh( X, Y, Z )
xs2pdf( gcf(), "images/10_mesh_v1.pdf" )

if getscilabmode() ~= "STD"
  quit()
end


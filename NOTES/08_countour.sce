function z = fun1(x, y)
  z = sin(x)*cos(y)
endfunction

function z = fun2(x, y)
  z = x^2 + y^2
endfunction

x = linspace(-%pi, %pi, 20)
y = linspace(-%pi, %pi, 20)
z1 = feval( x, y, fun1 )
z2 = feval( y, y, fun2 )

disp(size(z1))

clf()
N_level = 10
contour( x, y, z1, N_level ) // Alternatively: contour( x, y, my_surface, 10)
title("A contour plot")
xlabel('x')
ylabel('y')
xs2pdf( gcf(), "images/08_contour_v1.pdf" )


clf()
N_level = 10
contour( x, y, z2, N_level ) // Alternatively: contour( x, y, my_surface, 10)
title("A contour plot")
xlabel('x')
ylabel('y')
xs2pdf( gcf(), "images/08_contour_v2.pdf" )

clf()
plot3d(x, y, z1)
contour( x, y, z1, N_level, flag=[0 2 4])
xs2pdf( gcf(), "images/08_contour_v3.pdf" )


if getscilabmode() ~= "STD"
  quit()
end

x = linspace(-1, 1, 20)
y = linspace(-1, 1, 20)
[X,Y] = meshgrid(x,y)
fx = 0.5 .* X'
fy = 3 .* Y'

clf()
champ( x, y, fx, fy )
title("Vector field plot using champ")
xlabel('x')
ylabel('y')
xs2pdf( gcf(), "images/07_champ_v1.pdf" )

clf()
champ1( x, y, fx, fy )
title("Vector field plot using champ1")
xlabel('x')
ylabel('y')
xs2pdf( gcf(), "images/07_champ_v2.pdf" )

if getscilabmode() ~= "STD"
  quit()
end



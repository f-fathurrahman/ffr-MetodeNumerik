x = linspace( -10, 10, 50 )
y1 = x.^3 - 1.2*x.^2 - x

plot( x, y1 )
xs2pdf( gcf(), "images/01_plot_v1.pdf" )

plot( x, y1, "r*")
xs2pdf( gcf(), "images/01_plot_v2.pdf" )

L = 5
y2 = 100*cos(2*%pi/L*x)
plot( x, y1, 'r', x, y2)
xs2pdf( gcf(), "images/01_plot_v3.pdf" )


if getscilabmode() ~= "STD"
  quit()
end
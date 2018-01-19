x = linspace(-1,1,50)
T = 1.0
y = sin(2*%pi/T*x)

clf()
plot2d2( x, y )
xs2pdf( gcf(), "images/02_plot_v1.pdf" )

clf()
plot2d3( x, y )
xs2pdf( gcf(), "images/02_plot_v2.pdf" )

clf()
plot2d4( x, y )
xs2pdf( gcf(), "images/02_plot_v3.pdf" )

if getscilabmode() ~= "STD"
  quit()
end

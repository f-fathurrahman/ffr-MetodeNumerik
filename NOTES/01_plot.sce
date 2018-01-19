x = linspace( -3, 3, 50 )
y1 = x.^3 - 1.2*x.^2 - x

clf()
plot( x, y1 )
xs2pdf( gcf(), "images/01_plot_v1.pdf" )

clf()
plot( x, y1, "r*")
xs2pdf( gcf(), "images/01_plot_v2.pdf" )

L = 5
A = 10
y2 = 10*cos(2*%pi/L*x)
clf()
plot( x, y1, "r*", x, y2, "r-")
xtitle("This is a title")
xlabel("My x-label")
ylabel("My y-label")
xs2pdf( gcf(), "images/01_plot_v3.pdf" )


if getscilabmode() ~= "STD"
  quit()
end

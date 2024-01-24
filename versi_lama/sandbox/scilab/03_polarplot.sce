theta = 0:0.02:3*%pi

clf()
R = 1.2
polarplot( sin(3*theta), cos(2*theta))
xs2pdf( gcf(), "images/03_polar_plot_v1.pdf" )

if getscilabmode() ~= "STD"
  quit()
end

x = [1 2 3 4; 5 4 3 6]

clf()
Matplot(x)

xs2pdf( gcf(), "images/05_matplot_v1.pdf" )

if getscilabmode() ~= "STD"
  quit()
end

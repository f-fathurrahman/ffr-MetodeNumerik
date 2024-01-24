
mydata = rand(1,10e4,"normal")

N_bins = 30
histplot( N_bins, mydata )
title("An example of histogram of normally distributed data")
xlabel("Value")
ylabel("Probability")

xs2pdf( gcf(), "images/04_histogram_v1.pdf" )

if getscilabmode() ~= "STD"
  quit()
end

function PlotData(data,input)

imagesc(1:input.NR,1:input.NS, (data)); 
    xlabel('Receiver number p \rightarrow'); 
    ylabel('Source number q \rightarrow');  
    axis('equal','tight'); axis xy;
    colorbar('hor'); colormap jet;
function IMAGESC(x1,x2,matrix2D)

% interchange of x1 and x2
  imagesc(x2,x1,matrix2D);
  
% add some standard options
  xlabel('x_2 \rightarrow'); 
  ylabel('\leftarrow x_1');  
  axis('equal','tight');
  colorbar('hor'); colormap jet;



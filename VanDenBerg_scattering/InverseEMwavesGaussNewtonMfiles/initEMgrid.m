function input = initEMgrid(input)

   input.N1 = 50;                              % number of samples in x_1  
   input.N2 = 50;                              % number of samples in x_2 
   input.dx = 3;                                % with meshsize dx
   x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx;   
   x2 = -(input.N2+1)*input.dx/2 + (1:input.N2)*input.dx;
   [input.X1,input.X2] = ndgrid(x1,x2);
 
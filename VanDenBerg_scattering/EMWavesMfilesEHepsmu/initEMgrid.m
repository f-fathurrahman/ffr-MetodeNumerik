function input = initEMgrid(input)
global nDIM;  

if nDIM == 2;   % Grid in two-dimensional space ---------------------------

   input.N1 = 240;                              % number of samples in x_1  
   input.N2 = 200;                              % number of samples in x_2 
   input.dx = 1;                                % with meshsize dx
   x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx;   
   x2 = -(input.N2+1)*input.dx/2 + (1:input.N2)*input.dx;
   [input.X1,input.X2] = ndgrid(x1,x2);
 
elseif nDIM == 3;  % Grid in three-dimensional space ----------------------

   input.N1 = 240;                              % number of samples in x_1  
   input.N2 = 200;                              % number of samples in x_2
   input.N3 = 161;                              % number of samples in x_3
   input.dx = 1;                                % with meshsize dx
   x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx;   
   x2 = -(input.N2+1)*input.dx/2 + (1:input.N2)*input.dx;
   x3 = -(input.N3+1)*input.dx/2 + (1:input.N3)*input.dx;
   [input.X1,input.X2,input.X3] = ndgrid(x1,x2,x3);

end %if

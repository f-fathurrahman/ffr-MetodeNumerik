function input = initGrid(input)
global nDIM;  

if nDIM == 1      % Grid in one-dimensional space ----------------------- 
         
   input.N1 = 120;                              % number of samples in x_1 
   input.dx = 2;                                % meshsize dx 
   x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx; 
   input.X1 = x1';
   % Now X1 is a column vector equivalent with x1 axis pointing downwards    

elseif nDIM == 2   % Grid in two-dimensional space -----------------------

   input.N1 = 120;                              % number of samples in x_1  
   input.N2 = 100;                              % number of samples in x_2 
   input.dx = 2;                                % with meshsize dx
   x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx;   
   x2 = -(input.N2+1)*input.dx/2 + (1:input.N2)*input.dx;
   [input.X1,input.X2] = ndgrid(x1,x2);
   % Now array subscripts are equivalent with Cartesian coordinates
   % x1 axis points downwards and x2 axis is in horizontal direction
   % x1 = X1(:,1) is a column vector in vertical direction
   % x2 = X2(1,:) is a row vector in horizontal direction
 
elseif nDIM == 3   % Grid in three-dimensional space ----------------------

   input.N1 = 120;                              % number of samples in x_1  
   input.N2 = 100;                              % number of samples in x_2
   input.N3 = 80;                               % number of samples in x_3
   input.dx = 2;                                % with meshsize dx
   x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx;   
   x2 = -(input.N2+1)*input.dx/2 + (1:input.N2)*input.dx;
   x3 = -(input.N3+1)*input.dx/2 + (1:input.N3)*input.dx;
   [input.X1,input.X2,input.X3] = ndgrid(x1,x2,x3);
   % Now array subscripts are equivalent with Cartesian coordinates
   % x1 axis points downwards and x2 axis is in horizontal direction
   % x1 = X1(:,1,1) is a column vector in vertical direction
   % x2 = X2(1,:,1) is a row vector in first horizontal direction
   % x3 = X3(1,1,:) is a row vector in second horizontal direction
   
end %if

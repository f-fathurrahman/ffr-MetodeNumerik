function [input] = initContrast(input)
input.a  = 40;               % radius circle cylinder 
contrast = 1 - input.c_0^2/input.c_sct^2; 

if contrast > 0;  input.signCHI =  1; end
if contrast < 0;  input.signCHI = -1; end
input.contrast = contrast;

input.xO(1) = input.a / 2;   % center coordinates of circle
input.xO(2) = input.a / 3; 
          R = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2);          
  input.CHI = contrast .* (R < input.a);   
      
% Figures of contrast
x1 = input.X1(:,1);   x2 = input.X2(1,:); 
figure(1);  IMAGESC(x1,x2,real(input.CHI));  
figure(2);  mesh(abs(input.CHI),'edgecolor', 'k'); 
            view(37.5,45); axis xy; axis('off'); axis('tight') 
            axis([1 input.N1 1 input.N2  0 1.5*max(abs(input.CHI(:)))])
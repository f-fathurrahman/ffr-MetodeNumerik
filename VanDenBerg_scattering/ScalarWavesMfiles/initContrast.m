function [input] = initContrast(input)
global nDIM; 

input.a  = 40;   % half width slab / radius circle cylinder / radius sphere

contrast = 1 - input.c_0^2/input.c_sct^2;  

if     nDIM == 1;   R = sqrt(input.X1.^2);     
elseif nDIM == 2;   R = sqrt(input.X1.^2 + input.X2.^2);
elseif nDIM == 3;   R = sqrt(input.X1.^2 + input.X2.^2 + input.X3.^2);
end % if

input.CHI = contrast .* (R < input.a); 

function [Pn,Pn_1,Pn_2] = Legendre(n,Z,Pn,Pn_1,Pn_2)
% Computation of Legendre polynomials by recursion 
   if n == 0; Pn  = ones(size(Z)); Pn_2 = Pn; end;
   if n == 1; Pn  = Z;   end; 
   if n > 1;                   
       Pn_1 = Pn;              
       Pn = (2*n-1)/n * Z .* Pn_1 -(n-1)/n * Pn_2;  
       Pn_2 = Pn_1;
   end
end


function [adjGMdw] = AdjDOPMdw(f,input)
gam0 = input.gamma_0;    dx = input.dx;     xR = input.xR;
X1   = input.X1;         X2 = input.X2; 

adjGMdw    = cell(1,2); 
adjGMdw{1} = zeros(input.N1,input.N2);  
adjGMdw{2} = zeros(input.N1,input.N2);
for p = 1: input.NR   
    DIS      = sqrt((xR(1,p)-(X1+dx/2)).^2 + (xR(2,p)-X2).^2);
    d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d1_G     = ((xR(1,p)-(X1+dx/2))./DIS) .* d_G; 
  adjGMdw{1} = adjGMdw{1} + 2 * dx * conj(d1_G) * f(p); 
  
    DIS      = sqrt((xR(1,p)-X1).^2 + (xR(2,p)-(X2+dx/2)).^2);
    d_G      = - gam0 * (1/(2*pi)) .* besselk(1,gam0*DIS);
    d2_G     = ((xR(2,p)-(X2+dx/2))./DIS) .* d_G; 
  adjGMdw{2} = adjGMdw{2} + 2 * dx * conj(d2_G) * f(p);
end % p_loop                                             
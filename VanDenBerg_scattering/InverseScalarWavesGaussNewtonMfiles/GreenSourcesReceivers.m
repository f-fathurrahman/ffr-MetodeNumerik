function [G_S,G_R] = GreenSourcesReceivers(input) 

gam0 = input.gamma_0;   
X1   = input.X1; 
X2   = input.X2;
dx   = input.dx;
  
delta  = (pi)^(-1/2) * dx;           % radius circle with area of dx^2  
factor = 2 * besseli(1,gam0*delta) / (gam0*delta);
                                     % factor for weak form if DIS > delta 
 xS   = input.xS;
 G_S  = cell(1,input.NS);
for q = 1 : input.NS
     DIS   = sqrt( (X1-xS(1,q)).^2 + (X2-xS(2,q)).^2 );
    G_S{q} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);          
end %q_loop

  xR  = input.xR;
  G_R = cell(1,input.NR);
for p = 1 : input.NR    
     DIS  = sqrt((xR(1,p)-X1).^2 +(xR(2,p)-X2).^2);
   G_R{p} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);      
end % p_loop
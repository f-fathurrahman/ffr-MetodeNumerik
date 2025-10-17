function [G_S,dG_S,G_R,dG_R] = dGreenSourcesReceivers(input) 

gam0 = input.gamma_0;   X1 = input.X1;  X2 = input.X2;
dx   = input.dx;
  
delta  = (pi)^(-1/2) * dx;           % radius circle with area of dx^2  
factor = 2 * besseli(1,gam0*delta) / (gam0*delta);
                                     % factor for weak form if DIS > delta                                   
xS  = input.xS;
G_S = cell(1,input.NS); dG_S  = cell(2,input.NS);
for q = 1 : input.NS
      DIS   = sqrt((X1-xS(1,q)).^2 + (X2-xS(2,q)).^2);
     G_S{q} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);  
        dG  = - gam0 * factor * 1/(2*pi).* besselk(1,gam0*DIS);
  dG_S{1,q} = (X1-xS(1,q))./DIS .* dG;
  dG_S{2,q} = (X2-xS(2,q))./DIS .* dG; 
end %q_loop

xR  = input.xR;
G_R = cell(1,input.NR); dG_R = cell(2,input.NR);
for p = 1 : input.NR    
       DIS  = sqrt((xR(1,p)-X1).^2 +(xR(2,p)-X2).^2);
     G_R{p} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);  
        dG  = - gam0 * factor * 1/(2*pi).* besselk(1,gam0*DIS);
  dG_R{1,p} = - (xR(1,p)-X1)./DIS .* dG;
  dG_R{2,p} = - (xR(2,p)-X2)./DIS .* dG; 
end % p_loop
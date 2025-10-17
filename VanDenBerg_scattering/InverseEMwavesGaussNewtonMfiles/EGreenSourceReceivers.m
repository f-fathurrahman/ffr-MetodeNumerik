function pG = EGreenSourceReceivers(input) 

gam0 = input.gamma_0;   
dx   = input.dx;
delta  = (pi)^(-1/2) * dx;           % radius circle with area of dx^2  
factor = 2 * besseli(1,gam0*delta) / (gam0*delta);
                                     % factor for weak form if DIS > delta    
% Sources
 xS  = input.xS;
   G = cell(1,input.NS); 
dG11 = cell(1,input.NS); dG21 = cell(1,input.NS); dG22 = cell(1,input.NS);
for q = 1 : input.NS
      X1  = input.X1 - xS(1,q); X2 = input.X2 - xS(2,q);   
      DIS = sqrt(X1.^2 + X2.^2);  
       X1 = X1./DIS;  X2 = X2./DIS; 
     G{q} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);  
       dG = - factor * gam0 .* 1/(2*pi).* besselk(1,gam0*DIS);
  dG11{q} = (2*X1.*X1 - 1) .* (-dG./DIS) + gam0^2 * X1.*X1 .* G{q};
  dG21{q} =  2*X2.*X1      .* (-dG./DIS) + gam0^2 * X2.*X1 .* G{q};  
  dG22{q} = (2*X2.*X2 - 1) .* (-dG./DIS) + gam0^2 * X2.*X2 .* G{q};
end %q_loop
pG.GS = G;   pG.dGS11 = dG11;   pG.dGS21 = dG21;   pG.dGS22 = dG22; 

% Receivers
 xR  = input.xR;
   G = cell(1,input.NR); 
dG11 = cell(1,input.NR); dG21 = cell(1,input.NR); dG22 = cell(1,input.NR);
for p = 1 : input.NR    
      X1  = input.X1 - xR(1,p);  X2 = input.X2 - xR(2,p);   
      DIS = sqrt(X1.^2 + X2.^2); 
       X1 = X1./DIS;  X2 = X2./DIS;
     G{p} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);  
       dG = - factor * gam0 .* 1/(2*pi).* besselk(1,gam0*DIS);
  dG11{p} = (2*X1.*X1 - 1) .* (-dG./DIS) + gam0^2 * X1.*X1 .* G{p};
  dG21{p} =  2*X2.*X1      .* (-dG./DIS) + gam0^2 * X2.*X1 .* G{p};
  dG22{p} = (2*X2.*X2 - 1) .* (-dG./DIS) + gam0^2 * X2.*X2 .* G{p};
end % p_loop
pG.GR = G;   pG.dGR11 = dG11;   pG.dGR21 = dG21;   pG.dGR22 = dG22; 
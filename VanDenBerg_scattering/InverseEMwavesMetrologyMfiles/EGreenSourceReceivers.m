function pG = EGreenSourceReceivers(input) 

gam0  = input.gamma_0;      dx = input.dx; 
delta = (pi)^(-1/2) * dx;   % radius circle with area of dx^2  
factor = 2 * besseli(1,gam0*delta) / (gam0*delta);  % factor for weak form 
                                     
 xS  = input.xS;   % Sources
   G = cell(1,input.NS); dG1  = cell(1,input.NS); dG2  = cell(1,input.NS);
dG11 = cell(1,input.NS); dG21 = cell(1,input.NS); dG22 = cell(1,input.NS);

for q = 1 : input.NS
      X1  = input.X1 - xS(1,q); X2 = input.X2 - xS(2,q);   
      DIS = sqrt(X1.^2 + X2.^2);  
      X1  = X1./DIS;  X2 = X2./DIS; 
     G{q} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);  
     dG   = - factor * gam0 .* 1/(2*pi).* besselk(1,gam0*DIS);
   dG1{q} = X1 .* dG;
   dG2{q} = X2 .* dG;
  dG11{q} = (2*X1.*X1 - 1) .* (-dG./DIS) + gam0^2 * X1.*X1 .* G{q};
  dG21{q} =  2*X2.*X1      .* (-dG./DIS) + gam0^2 * X2.*X1 .* G{q};  
  dG22{q} = (2*X2.*X2 - 1) .* (-dG./DIS) + gam0^2 * X2.*X2 .* G{q};
end %q_loop
pG.GS    = G;      pG.dGS1  = dG1;    pG.dGS2  = dG2;   
pG.dGS11 = dG11;   pG.dGS21 = dG21;   pG.dGS22 = dG22; 

 xR  = input.xR; % Receivers
   G = cell(1,input.NR); dG1  = cell(1,input.NR); dG2  = cell(1,input.NR);
dG11 = cell(1,input.NR); dG21 = cell(1,input.NR); dG22 = cell(1,input.NR);
for p = 1 : input.NR    
      X1  = input.X1 - xR(1,p);  X2 = input.X2 - xR(2,p);   
      DIS = sqrt(X1.^2 + X2.^2); 
       X1 = X1./DIS;  X2 = X2./DIS;
     G{p} = factor * 1/(2*pi) .* besselk(0,gam0*DIS);  
     dG   = - factor * gam0 .* 1/(2*pi).* besselk(1,gam0*DIS);
   dG1{p} = X1 .* dG;
   dG2{p} = X2 .* dG;    
  dG11{p} = (2*X1.*X1 - 1) .* (-dG./DIS) + gam0^2 * X1.*X1 .* G{p};
  dG21{p} =  2*X2.*X1      .* (-dG./DIS) + gam0^2 * X2.*X1 .* G{p};
  dG22{p} = (2*X2.*X2 - 1) .* (-dG./DIS) + gam0^2 * X2.*X2 .* G{p};
end; % p_loop
pG.GR    = G;      pG.dGR1  = dG1;    pG.dGR2  = dG2;  
pG.dGR11 = dG11;   pG.dGR21 = dG21;   pG.dGR22 = dG22;
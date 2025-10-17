function [E_inc] = IncEMwave(pG,q,input) 

gam0 = input.gamma_0;    xS = input.xS;  
  M1 = xS(1,q) / sqrt(xS(1,q)^2 + xS(2,q)^2);   
  M2 = xS(2,q) / sqrt(xS(1,q)^2 + xS(2,q)^2);
      
E_inc{1} = (-gam0^2 * pG.GS{q} + pG.dGS11{q}) * M1 + pG.dGS21{q} * M2;
E_inc{2} = (-gam0^2 * pG.GS{q} + pG.dGS22{q}) * M2 + pG.dGS21{q} * M1;
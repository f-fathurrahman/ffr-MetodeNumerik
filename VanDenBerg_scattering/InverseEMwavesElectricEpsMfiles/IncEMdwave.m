function [E_inc,dEinc] = IncEMdwave(pG,q,input) 
gam0 = input.gamma_0;    xS = input.xS;  

% Non shifted grid (0) --------------------------------------------------
        M1 = xS(1,q) / sqrt(xS(1,q)^2 + xS(2,q)^2);   
        M2 = xS(2,q) / sqrt(xS(1,q)^2 + xS(2,q)^2);
  E_inc{1} = (-gam0^2*pG.GS{q} + pG.dGS11{q}) * M1 + pG.dGS21{q} * M2;
  E_inc{2} = (-gam0^2*pG.GS{q} + pG.dGS22{q}) * M2 + pG.dGS21{q} * M1;
   
% Shifted grid (1) -----------------------------------------------------
  dEinc{1} = 0 * E_inc{1};   
         i = 1:input.N1-1;
  dEinc{1}(i,:) = (E_inc{1}(i+1,:) + E_inc{1}(i,:))/2;
   
% Shifted grid (2) ------------------------------------------------------   
  dEinc{2} = 0 * E_inc{2};   
         j = 1:input.N2-1;
  dEinc{2}(:,j) = (E_inc{2}(:,j+1) + E_inc{2}(:,j))/2;
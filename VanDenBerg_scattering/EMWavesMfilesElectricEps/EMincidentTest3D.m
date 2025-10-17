clear all; clc; close all; clear workspace
input = initEM();  gam0 = input.gamma_0;  xS = input.xS;
 
% (1) Cartesian coordinates: field from electric dipole in negative x_1 ---
   X1   = input.X1-xS(1);   X2 = input.X2-xS(2);    X3 = input.X3-xS(3); 
  DIS   = sqrt(X1.^2 + X2.^2 + X3.^2);
   X1   = X1 ./ DIS;        X2 = X2 ./ DIS;         X3 = X3 ./ DIS;
    G   = exp(-gam0*DIS) ./ (4*pi*DIS);
    dG  = -(gam0 + 1./DIS) .* G;
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2 *X1.*X1) .* G;
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2 *X2.*X1) .* G;
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2 *X3.*X1) .* G;
   
  E{1}  = - (-gam0^2 * G + dG11);
  E{2}  = - dG21;
  E{3}  = - dG31;
  
  % H1 is zero !
  ZH{2} = -gam0 * X3 .* dG;
  ZH{3} =  gam0 * X2 .* dG;

% (2) Compute the field in terms  of Bessel series ------------------------
  [E1,E2,E3,ZH2,ZH3] = IncEMsphere(input);
  
% (3) Plot symmetric fields with respect to x3=0 at (x3=0 or x3=dx/2) -----   
  N3cross = floor(input.N3/2+1); 
  x1 = input.X1(:,1,N3cross);  x2 = input.X2(1,:,N3cross);  

  plotFieldError('3D: abs(E_1)', E{1}, E1, x1,x2,input);
  plotFieldError('3D: abs(E_2)', E{2}, E2, x1,x2,input);
  plotFieldError('3D: abs(Z_0H_3)', ZH{3}, ZH3, x1,x2,input); 
clear all; clc; close all; clear workspace
input = initEM();  gam0 = input.gamma_0;      
    
% (1) Cartesian coordinates: electric dipole in negative x_1 --------------

  xS   = input.xS;
  X1   = input.X1-xS(1);     
  X2   = input.X2-xS(2); 
  DIS  = sqrt(X1.^2 + X2.^2); X1 = X1./DIS;  X2 = X2./DIS;
   G   = 1/(2*pi).* besselk(0,gam0*DIS);  
  dG   = - gam0 .* 1/(2*pi).* besselk(1,gam0*DIS);
  dG11 = (2*X1.*X1 - 1) .* (-dG./DIS) + gam0^2 * X1.*X1 .* G;
  dG21 =  2*X2.*X1      .* (-dG./DIS) + gam0^2 * X2.*X1 .* G;
  
  E{1} = - (-gam0^2 * G + dG11);
  E{2} = - dG21; 
 ZH{3} = gam0 * X2 .* dG; 

% (2) Transform Cartesian coordinates to polar ccordinates ----------------

  rS = sqrt(xS(1)^2+xS(2)^2);      phiS = atan2(xS(2),xS(1));  
  X1 = input.X1;   
  X2 = input.X2;
  R  = sqrt(X1.^2+X2.^2+1e-16);    PHI  = atan2(X2,X1);  
  
  % and compute Bessel series with M terms --------------------------------
  M = 50;                                   % increase M for more accuracy
  Er = zeros(size(R));  Ephi = zeros(size(R));  ZH3 = zeros(size(R));  
  for m = 1 : M
    arg0 = gam0*R;       Ib0 = besseli(m,arg0);    
                        dIb0 = besseli(m+1,arg0) + m./arg0 .* Ib0;  
                         KbS = besselk(m,gam0*rS);
     Er  = Er   + 2*m^2 .* Ib0 .*  KbS .* cos(m*(PHI-phiS)); 
    Ephi = Ephi - 2*m   .*dIb0 .*  KbS .* sin(m*(PHI-phiS));
    ZH3  = ZH3  - 2*m   .* Ib0 .*  KbS .* sin(m*(PHI-phiS));           
  end % m_loop
  
  Er   =         1/(2*pi) * Er./R ./rS;
  Ephi =  gam0 * 1/(2*pi) * Ephi  ./rS;
  ZH3  = -gam0 * 1/(2*pi) * ZH3   ./rS;
  
% (3) Determine mean error and plot error in Cartesian vectors ------------

  E1 = cos(PHI) .* Er - sin(PHI) .* Ephi;  
  E2 = sin(PHI) .* Er + cos(PHI) .* Ephi;
  
  x1 = input.X1(:,1);   x2 = input.X2(1,:);
  plotFieldError('2D: abs(E_1)',    E{1},  E1, x1,x2,input);         
  plotFieldError('2D: abs(E_2)',    E{2},  E2, x1,x2,input);  
  plotFieldError('2D: abs(ZH_3)',  ZH{3}, ZH3, x1,x2,input); 
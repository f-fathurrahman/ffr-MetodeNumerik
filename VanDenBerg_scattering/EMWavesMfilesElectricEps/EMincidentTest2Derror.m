clear all; clc; close all; clear workspace
input = initEM();  
gam0  = input.gamma_0;     
   a  = input.a;    
  PHI = 0:.01:2*pi; 
    
% (1) Cartesian coordinates: electric dipole in negative x_1 --------------

   xS  = input.xS; 
   x1  = a*cos(PHI)-xS(1);      x2 = a*sin(PHI)-xS(2); 
  DIS  = sqrt(x1.^2 + x2.^2);   x1 = x1./DIS;  x2 = x2./DIS;
   G   = 1/(2*pi).* besselk(0,gam0*DIS);  
  dG   = - gam0 .* 1/(2*pi).* besselk(1,gam0*DIS);
  dG11 = (2*x1.*x1 - 1) .* (-dG./DIS) + gam0^2 * x1.*x1 .* G;
  dG21 =  2*x2.*x1      .* (-dG./DIS) + gam0^2 * x2.*x1 .* G;
  E{1} = - (-gam0^2 * G + dG11);
  E{2} = - dG21; 
 ZH{3} = gam0 * x2 .* dG; 

% (2) Polar coordinates: expansion in bessel functions --------------------

  rS   = sqrt(xS(1)^2+xS(2)^2);   
  phiS = atan2(xS(2),xS(1));  
  Er   = zeros(size(PHI));  
  Ephi = zeros(size(PHI));  
  ZH3  = zeros(size(PHI));    
  
  M = 20;                                   % increase M for more accuracy
  for m = 1 : M  
    Ib0  = besseli(m,gam0*a);  dIb0 = besseli(m+1,gam0*a)+m*Ib0/(gam0*a);  
    KbS  = besselk(m,gam0*rS);
    Er   = Er   + 2*m^2 * Ib0 *  KbS .* cos(m*(PHI-phiS)); 
    Ephi = Ephi - 2*m   *dIb0 *  KbS .* sin(m*(PHI-phiS));
    ZH3  = ZH3  - 2*m   * Ib0 *  KbS .* sin(m*(PHI-phiS));           
  end % m_loop
  
  Er   =         1/(2*pi) * Er/a /rS;
  Ephi =  gam0 * 1/(2*pi) * Ephi /rS;
  ZH3  = -gam0 * 1/(2*pi) * ZH3  /rS;
  
  E1 = cos(PHI) .* Er - sin(PHI) .* Ephi;    
  E2 = sin(PHI) .* Er + cos(PHI) .* Ephi; 

% Print the normalized error on the circular boundary at r = a ----------  
  norm_Error = norm(E{1}(:)-E1(:),1) / norm(E{1}(:),1);      
               disp(['Error ',' E_1',' = ',num2str(norm_Error)]); 
  norm_Error = norm(E{2}(:)-E2(:),1) / norm(E{2}(:),1);      
               disp(['Error ',' E_2',' = ',num2str(norm_Error)]); 
  norm_Error = norm(ZH{3}(:)-ZH3(:),1) / norm(ZH{3}(:),1);      
               disp(['Error ','ZH_3',' = ',num2str(norm_Error)]); 
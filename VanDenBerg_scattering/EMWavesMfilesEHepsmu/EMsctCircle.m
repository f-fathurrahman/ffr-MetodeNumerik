clear all; clc; close all;
input = initEM();       
      a = input.a;   
   gam0 = input.gamma_0; 
eps_sct = input.eps_sct; 
 mu_sct = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct); 
 Z_sct  = sqrt(mu_sct/eps_sct);  

% (1) Transform Cartesian coordinates to polar coordinates ---------------- 
  xR = input.xR; 
  xS = input.xS; 
  rR = sqrt(xR(1,:).^2 + xR(2,:).^2);  phiR = atan2(xR(2,:),xR(1,:));  
  rS = sqrt(xS(1)^2+xS(2)^2);          phiS = atan2(xS(2),xS(1)); 

% (2) Compute coefficients of Bessel series expansion ---------------------
  arg0 = gam0 * input.a;  args = gam_sct * input.a; 
  M = 20;    A = zeros(1,M);                % increase M for more accuracy 
  for m = 1 : M;                      
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0;
    denominator = Z_sct * dIbs*Kb0 - dKb0*Ibs;
    A(m)  =     -(Z_sct * dIbs*Ib0 - dIb0*Ibs) / denominator; 
  end
  
% (3) Compute reflected Er field at receivers (data) ----------------------   
  Er = zeros(size(rR)); Ephi = zeros(size(rR));  ZH3 = zeros(size(rR));
  for m = 1 : M;
    arg0 = gam0*rR;  Kb0 =  besselk(m,arg0);   
                    dKb0 = -besselk(m+1,arg0) + m./arg0 .* Kb0; 
                     KbS = besselk(m,gam0*rS);
    Er   = Er   + A(m)*2*m^2.* Kb0.* KbS .*cos(m*(phiR-phiS));
    Ephi = Ephi - A(m)*2*m  .*dKb0.* KbS .*sin(m*(phiR-phiS));
    ZH3  = ZH3  - A(m)*2*m  .* Kb0.* KbS .*sin(m*(phiR-phiS)); 
  end % m_loop

  Er   =         1/(2*pi) * Er./rR ./rS;
  Ephi =  gam0 * 1/(2*pi) * Ephi   ./rS;
  ZH3  = -gam0 * 1/(2*pi) * ZH3    ./rS;
  E{1} =   cos(phiR) .* Er  - sin(phiR) .* Ephi;
  E{2} =   sin(phiR) .* Er  + cos(phiR) .* Ephi;
  Edata2D = sqrt(abs(E{1}).^2 + abs(E{2}).^2); 
  Hdata2D = abs(ZH3);
  
  if exist(fullfile(cd, 'EDATA2D.mat'), 'file'); delete EDATA2D.mat; end;
  if exist(fullfile(cd, 'HDATA2D.mat'), 'file'); delete HDATA2D.mat; end;
  
  displayEdata(Edata2D,input);                   save EDATA2D Edata2D;
  displayHdata(Hdata2D,input);                   save HDATA2D Hdata2D;
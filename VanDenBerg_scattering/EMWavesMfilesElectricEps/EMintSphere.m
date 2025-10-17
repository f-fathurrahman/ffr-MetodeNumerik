clear all; clc; close all; clear workspace
input = initEM();   a = input.a;   gam0 = input.gamma_0;  xS = input.xS;
eps_sct = input.eps_sct;               mu_sct = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct); Z_sct  = sqrt(mu_sct/eps_sct); 

% (1) Transform Cartesian coordinates to polar ccordinates ---------------- 
   rS   = sqrt(xS(1)^2+xS(2)^2+xS(3)^2);     
  phiS = atan2(xS(2),xS(1));            thetaS = acos(xS(3)/rS);
  X1   = input.X1;       X2 = input.X2;     X3 = input.X3;
  R    = sqrt(X1.^2+X2.^2+X3.^2+1e-16);    
  PHI  = atan2(X2,X1);                  THETA  = acos(X3./R);
% (2) Compute coefficients of series expansion ----------------------------
  N = 20;                                 % increase N for more accuracy 
  A = zeros(1,N);  B = zeros(1,N);
  arg0 = gam0 * input.a;  args = gam_sct * input.a;
  for n = 1 : N
    Ib0 = besseli(n+1/2,arg0);  dIb0 = besseli(n+3/2,arg0)+(n+1)/arg0*Ib0; 
    Ibs = besseli(n+1/2,args);  dIbs = besseli(n+3/2,args)+(n+1)/args*Ibs; 
    Kb0 = besselk(n+1/2,arg0);  dKb0 =-besselk(n+3/2,arg0)+(n+1)/arg0*Kb0;
    A(n) = -(Z_sct *dIbs*Ib0 - dIb0*Ibs) / (Z_sct *dIbs*Kb0 - dKb0*Ibs);
  % factor sqrt(pi^2/4/argo/args) in numerator and denominator is omitted
    B(n) = sqrt(args/arg0) * (Ib0 + A(n)*Kb0) / Ibs;  
  end
% (3) Compute interior wave field using N terms of Bessel series expansion
  COS     =  cos(THETA)*cos(thetaS)+sin(THETA)*sin(thetaS).*cos(PHI-phiS);
  dTHETA  = -sin(THETA)*cos(thetaS)+cos(THETA)*sin(thetaS).*cos(PHI-phiS);  
  dPHI_sin=                                  - sin(thetaS).*sin(PHI-phiS);
  Pn  = zeros(size(R));  Pn_1 = Pn;  Pn_2 = Pn;    dPn  = Pn;
  Er  = zeros(size(R));  Etheta = zeros(size(R));  Ephi = zeros(size(R));
  [Pn,Pn_1,Pn_2] = Legendre(0,COS,Pn,Pn_1,Pn_2);
  for n = 1 : N  
    dPn = n*Pn + COS.*dPn; [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);                   
    arg = gam_sct*R; fctr =sqrt(pi/2./arg); Ibs = fctr.*besseli(n+1/2,arg); 
                     dIbs = fctr .* besseli(n+3/2,arg) + (n+1)./arg.*Ibs;
    arg = gam0*rS;   fctr =sqrt(pi/2 /arg); KbS = fctr *besselk(n+1/2,arg);
    Er     = Er     + B(n) *n*(n+1)*(2*n+1).* Ibs .* KbS            .* Pn; 
    Etheta = Etheta + B(n) *        (2*n+1).*dIbs .* KbS .* dTHETA  .*dPn;
    Ephi   = Ephi   + B(n) *        (2*n+1).*dIbs .* KbS .* dPHI_sin.*dPn;                   
  end % n_loop 
  Er     = (gam0        /eps_sct) * Er     ./ (2*pi^2 .*R*rS);
  Etheta = (gam0*gam_sct/eps_sct) * Etheta ./ (2*pi^2    *rS);
  Ephi   = (gam0*gam_sct/eps_sct) * Ephi   ./ (2*pi^2    *rS);
% (4) Transform to Cartesian vectors --------------------------------------
E1 = sin(THETA).*cos(PHI).*Er+cos(THETA).*cos(PHI).*Etheta-sin(PHI).*Ephi;
E2 = sin(THETA).*sin(PHI).*Er+cos(THETA).*sin(PHI).*Etheta+cos(PHI).*Ephi;
E3 = cos(THETA)          .*Er-sin(THETA)          .*Etheta               ;
w_E{1} = input.CHI_eps .* E1;
w_E{2} = input.CHI_eps .* E2;
w_E{3} = input.CHI_eps .* E3;
plotContrastSourcewE(w_E,input);  save CONTRASTSOURCE w_E; 
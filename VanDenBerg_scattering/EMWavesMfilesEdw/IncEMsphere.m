function [E1,E2,E3,ZH2,ZH3] = IncEMsphere(input)
gam0 = input.gamma_0;  xS = input.xS;

% Transform Cartesian coordinates to spherical coordinates -----------
  rS   = sqrt(xS(1)^2+xS(2)^2+xS(3)^2);     
  phiS = atan2(xS(2),xS(1));            thetaS = acos(xS(3)/rS);
  X1   = input.X1;       X2 = input.X2;     X3 = input.X3;
  R    = sqrt(X1.^2+X2.^2+X3.^2+1e-16);    
  PHI  = atan2(X2,X1);                  THETA  = acos(X3./R);
% and compute incident wave as Bessel series with 0:N terms --------------
  N = 50;                                   % increase N for more accuracy
  COS    =  cos(THETA)*cos(thetaS)+sin(THETA)*sin(thetaS).*cos(PHI-phiS);
  dTHETA = -sin(THETA)*cos(thetaS)+cos(THETA)*sin(thetaS).*cos(PHI-phiS);
  dPHI_sin=                                 - sin(thetaS).*sin(PHI-phiS);     
  Pn  = zeros(size(R));   Pn_1 = Pn;  Pn_2 = Pn;   dPn = Pn;
  Er  = zeros(size(R));   Etheta = zeros(size(R));  Ephi = zeros(size(R)); 
                         ZHtheta = zeros(size(R)); ZHphi = zeros(size(R));
  [Pn,Pn_1,Pn_2] = Legendre(0,COS,Pn,Pn_1,Pn_2);  
  for n = 1 : N; 
    dPn = n*Pn + COS.*dPn; [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);
    arg = gam0*R;  fctr=sqrt(pi/2./arg); Ib0 = fctr.* besseli(n+1/2,arg);         
          dIb0 = fctr.*besseli(n+3/2,arg) + (n+1)./arg.*Ib0;
    arg = gam0*rS; fctr=sqrt(pi/2./arg); KbS = fctr * besselk(n+1/2,arg);
    Er     = Er + n*(n+1)*(2*n+1)* Ib0.* KbS             .*  Pn; 
    Etheta = Etheta      +(2*n+1)*dIb0.* KbS .* dTHETA   .* dPn;
    Ephi   = Ephi        +(2*n+1)*dIb0.* KbS .* dPHI_sin .* dPn;     
   ZHtheta = ZHtheta     +(2*n+1)* Ib0.* KbS .* dPHI_sin .* dPn;
   ZHphi   = ZHphi       -(2*n+1)* Ib0.* KbS .* dTHETA   .* dPn;
  end % n_loop 
   Er     = gam0   * Er      ./ (2*pi^2 .*R*rS);
   Etheta = gam0^2 * Etheta  ./ (2*pi^2    *rS); 
   Ephi   = gam0^2 * Ephi    ./ (2*pi^2    *rS); 
% ZHr     = 0   for TM case
  ZHtheta = gam0^2 * ZHtheta ./ (2*pi^2    *rS); 
  ZHphi   = gam0^2 * ZHphi   ./ (2*pi^2    *rS);
  
% Transform to Cartesian vectors
 E1 = sin(THETA).*cos(PHI).*Er+cos(THETA).*cos(PHI).*Etheta-sin(PHI).*Ephi;  
 E2 = sin(THETA).*sin(PHI).*Er+cos(THETA).*sin(PHI).*Etheta+cos(PHI).*Ephi;
 E3 = cos(THETA)          .*Er-sin(THETA)          .*Etheta;
 % ZH1= 0;
 ZH2 = cos(THETA).*sin(PHI).*ZHtheta + cos(PHI).*ZHphi;  
 ZH3 = -sin(THETA)         .*ZHtheta;
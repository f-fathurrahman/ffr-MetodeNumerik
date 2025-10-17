clear all; clc; close all; clear workspace; input = initEM(); a = input.a;
   gam0 = input.gamma_0;  eps_sct = input.eps_sct;  mu_sct = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct);      Z_sct  = sqrt(mu_sct/eps_sct); 
if exist(fullfile(cd, 'EDATA3D.mat'), 'file');   delete EDATA3D.mat;  end;
if exist(fullfile(cd, 'HDATA3D.mat'), 'file');   delete HDATA3D.mat;  end;
% (1) Transform Cartesian coordinates to polar coordinates ----------------
  xS   = input.xS;            rS = sqrt(xS(1)^2+xS(2)^2+xS(3)^2);     
  phiS = atan2(xS(2),xS(1));       thtS = acos(xS(3)/rS);
  xR   = input.xR;            rR = sqrt(xR(1,:).^2+xR(2,:).^2+xR(3,:).^2);
  phiR = atan2(xR(2,:),xR(1,:));   thtR = acos(xR(3,:)./rR(1,:));  
% (2) Compute Bessel series expansion at the reciever points --------------
  N = 20;                                 % increase N for more accuracy 
  A = zeros(1,N);  arg0 = gam0 * input.a;  args = gam_sct * input.a;
  for n = 1 : N;  
    Ib0 = besseli(n+1/2,arg0);  dIb0 = besseli(n+3/2,arg0)+(n+1)/arg0*Ib0; 
    Ibs = besseli(n+1/2,args);  dIbs = besseli(n+3/2,args)+(n+1)/args*Ibs; 
    Kb0 = besselk(n+1/2,arg0);  dKb0 =-besselk(n+3/2,arg0)+(n+1)/arg0*Kb0;
    A(n) = -(Z_sct *dIbs*Ib0 - dIb0*Ibs) / (Z_sct *dIbs*Kb0 - dKb0*Ibs);
  end  % factor sqrt(pi^2/4/argo/args) in numerator/denominator is omitted
  COS  = cos(thtR)*cos(thtS)+sin(thtR)*sin(thtS).*cos(phiR-phiS);
 dthtR =-sin(thtR)*cos(thtS)+cos(thtR)*sin(thtS).*cos(phiR-phiS);  
 dphiR_sin=                          - sin(thtS).*sin(phiR-phiS);     
  Pn  = zeros(size(rR));  Pn_1 = Pn;   Pn_2 = Pn;    dPn = Pn;
  Er  = zeros(size(rR));  Etht = Er;   Ephi = Er;  ZHtht = Er; ZHphi = Er;
  [Pn,Pn_1,Pn_2] = Legendre(0,COS,Pn,Pn_1,Pn_2);  
  for n = 1 : N; 
    dPn = n*Pn + COS.*dPn; [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);
    arg = gam0*rR; fctr = sqrt(pi/2./arg); Kb0 = fctr.*besselk(n+1/2,arg); 
                   dKb0 = - fctr .* besselk(n+3/2,arg) + (n+1)./arg.*Kb0;
    arg = gam0*rS; fctr = sqrt(pi/2 /arg); KbS = fctr *besselk(n+1/2,arg);
    Er   = Er    + A(n)*n*(n+1)*(2*n+1) * Kb0 .* KbS             .* Pn; 
    Etht = Etht  + A(n)        *(2*n+1) *dKb0 .* KbS .* dthtR    .*dPn;
    Ephi = Ephi  + A(n)        *(2*n+1) *dKb0 .* KbS .* dphiR_sin.*dPn;
   ZHtht = ZHtht + A(n)        *(2*n+1) * Kb0 .* KbS .* dphiR_sin.*dPn;
   ZHphi = ZHphi - A(n)        *(2*n+1) * Ib0 .* KbS .* dthtR    .*dPn;
  end % n_loop 
   Er   = gam0   * Er    ./ (2*pi^2 .*rR*rS);
   Etht = gam0^2 * Etht  ./ (2*pi^2     *rS);
   Ephi = gam0^2 * Ephi  ./ (2*pi^2     *rS);
  ZHtht = gam0^2 * ZHtht ./ (2*pi^2     *rS); 
  ZHphi = gam0^2 * ZHtht ./ (2*pi^2     *rS);
 E1 = sin(thtR).*cos(phiR).*Er+cos(thtR).*cos(phiR).*Etht-sin(phiR).*Ephi;
 E2 = sin(thtR).*sin(phiR).*Er+cos(thtR).*sin(phiR).*Etht+cos(phiR).*Ephi;
 E3 =            cos(thtR).*Er          - sin(thtR).*Etht;
ZH1 = cos(thtR).*cos(phiR).*ZHtht - sin(phiR).*ZHphi; 
ZH2 = cos(thtR).*sin(phiR).*ZHtht + cos(phiR).*ZHphi; 
ZH3 =           -sin(thtR).*ZHtht;
  Edata3D = sqrt(abs(E1).^2 + abs(E2).^2 + abs(E3).^2);
  Hdata3D = sqrt(abs(ZH1).^2 + abs(ZH2).^2 + abs(ZH3).^2);
%   Edata3D = E1+E2+E3;
%   Hdata3D = ZH1+ZH2+ZH3;

  displayEdata(Edata3D,input);             save EDATA3D Edata3D;
  displayHdata(Hdata3D,input);             save HDATA3D Hdata3D;
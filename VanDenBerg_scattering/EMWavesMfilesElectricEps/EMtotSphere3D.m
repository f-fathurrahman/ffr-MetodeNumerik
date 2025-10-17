%clear all; clc; close all; clear workspace
input = initEM();   a = input.a;   gam0 = input.gamma_0;  
eps_sct = input.eps_sct;               mu_sct = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct); Z_sct  = sqrt(mu_sct/eps_sct); 

% (1) Cartesian coordinates: field from electric dipole in negative x_1 --
 X1=input.X1-input.xS(1); X2=input.X2-input.xS(2); X3=input.X3-input.xS(3); 
  DIS   = sqrt(X1.^2 + X2.^2 + X3.^2+1e-16);
   X1   = X1 ./ DIS;        X2 = X2 ./ DIS;         X3 = X3 ./ DIS;
    G   = exp(-gam0*DIS) ./ (4*pi*DIS);  
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2 *X1.*X1) .* G;
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2 *X2.*X1) .* G;
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2 *X3.*X1) .* G;                            
  Einc1 = - (-gam0^2 * G + dG11);
  Einc2 = - dG21; 
  Einc3 = - dG31;
  
% Transform Cartesian coordinates to polar ccordinates --------------- 
    xS = input.xS;                    rS = sqrt(xS(1)^2+xS(2)^2+xS(3)^2);     
  phiS = atan2(xS(2),xS(1));      thetaS = acos(xS(3)/rS);
  X1   = input.X1;  X2 = input.X2;    X3 = input.X3;
  R    = sqrt(X1.^2+X2.^2+X3.^2+1e-16);    
  PHI  = atan2(X2,X1);            THETA  = acos(X3./R);

% (2) Compute coefficients of series expansion ---------------------------
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

% (3) Compute exterior, interior and total fields ------------------------- 
  COS     =  cos(THETA)*cos(thetaS)+sin(THETA)*sin(thetaS).*cos(PHI-phiS);
  dTHETA  = -sin(THETA)*cos(thetaS)+cos(THETA)*sin(thetaS).*cos(PHI-phiS);  
  dPHI_sin=                                  - sin(thetaS).*sin(PHI-phiS);
  Pn  = zeros(size(R));  Pn_1 = Pn;  Pn_2 = Pn;    dPn  = 0;
  Er_out=zeros(size(R)); Etheta_out=zeros(size(R)); Ephi_out=zeros(size(R));
  Er_int=zeros(size(R)); Etheta_int=zeros(size(R)); Ephi_int=zeros(size(R));

  [Pn,Pn_1,Pn_2] = Legendre(0,COS,Pn,Pn_1,Pn_2);  
  for n = 1 : N 
    dPn  = n*Pn + COS.*dPn; [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);
    arg = gam0*R;  fctr=sqrt(pi/2./arg); 
          Ib0 = fctr.* besseli(n+1/2,arg);        
          dIb0 = fctr.*besseli(n+3/2,arg) + (n+1)./arg.*Ib0;
           Kb0 = fctr.*besselk(n+1/2,arg);
          dKb0 = -fctr .* besselk(n+3/2,arg) + (n+1)./arg .*Kb0;
    arg = gam_sct*R; fctr =sqrt(pi/2./arg); Ibs = fctr.*besseli(n+1/2,arg); 
                     dIbs = fctr .* besseli(n+3/2,arg) + (n+1)./arg.*Ibs;

          
    arg = gam0*rS; fctr=sqrt(pi/2./arg); KbS = fctr * besselk(n+1/2,arg);
    
    
    Er_out     = Er_out + n*(n+1)*(2*n+1)* (Ib0+A(n)*Kb0) .* KbS            .*  Pn; 
    Etheta_out = Etheta_out      +(2*n+1)*(dIb0+A(n)*dKb0).* KbS .* dTHETA  .* dPn;
    Ephi_out   = Ephi_out        +(2*n+1)*(dIb0+A(n)*dKb0).* KbS .* dPHI_sin.* dPn; 
    
    Er_int     = Er_int + n*(n+1)*(2*n+1)* B(n).* Ibs.* KbS            .*  Pn; 
    Etheta_int = Etheta_int      +(2*n+1)* B(n).*dIbs.* KbS .* dTHETA  .* dPn;
    Ephi_int   = Ephi_int        +(2*n+1)* B(n).*dIbs.* KbS .* dPHI_sin.* dPn;  
  end % n_loop 
   Er_out(ceil(input.N1/2),ceil(input.N2/2),ceil(input.N3/2)) = 0; 
   Etheta_out(ceil(input.N1/2),ceil(input.N2/2),ceil(input.N3/2)) = 0;
   Ephi_out(ceil(input.N1/2),ceil(input.N2/2),ceil(input.N3/2)) = 0;
   
  Er    = (gam0                ) * Er_out    ./ (2*pi^2 .*R*rS) .* (R >= a-input.dx/2)...
         +(gam0        /eps_sct) * Er_int    ./ (2*pi^2 .*R*rS) .* (R <  a-input.dx/2);
  Etheta= (gam0*gam0           ) * Etheta_out./ (2*pi^2    *rS) .* (R >= a-input.dx/2)...
         +(gam0*gam_sct/eps_sct) * Etheta_int./ (2*pi^2    *rS) .* (R <  a-input.dx/2);
  Ephi  = (gam0*gam0           ) * Ephi_out  ./ (2*pi^2    *rS) .* (R >= a-input.dx/2)...
         +(gam0*gam_sct/eps_sct) * Ephi_int  ./ (2*pi^2    *rS) .* (R <  a-input.dx/2);
  % Transform to Cartesian vectors -----------------------------------------
E{1} = sin(THETA).*cos(PHI).*Er+cos(THETA).*cos(PHI).*Etheta-sin(PHI).*Ephi;
E{2} = sin(THETA).*sin(PHI).*Er+cos(THETA).*sin(PHI).*Etheta+cos(PHI).*Ephi;
E{3} = cos(THETA)          .*Er-sin(THETA)          .*Etheta               ;
  
% Scattered field
 E{1} =E{1} - Einc1;  
 E{2} =E{2} - Einc2; 
 E{3} =E{3} - Einc3;


    N1 = input.N1;         N2 = input.N2;    N3cross = floor(input.N3/2+1);
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1);   
   phi = 0:.01:2*pi; a =input.a
set(figure,'Units','centimeters','Position',[5 5 18 12]);             
subplot(1,3,1)     
IMAGESC(x1,x2,abs(E{1}(3:N1-2,3:N2-2,N3cross)));  hold on; 
    title('\fontsize{13} abs(E_1^{sct})');  plot(a*cos(phi),a*sin(phi),'w'); 
subplot(1,3,2)   
IMAGESC(x1,x2,abs(E{2}(3:N1-2,3:N2-2,N3cross)));  hold on; 
    title('\fontsize{13} abs(E_2^{sct})');  plot(a*cos(phi),a*sin(phi),'w');
subplot(1,3,3)
IMAGESC(x1,x2,abs(E{3}(3:N1-2,3:N2-2,N3cross)));   hold on; 
title('\fontsize{13} abs(E_3^{sct})');   plot(a*cos(phi),a*sin(phi),'w'); 

save SCATTEREDFIELD  E;   
 

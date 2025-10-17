clear all; clc; close all;
input = init(); 
 c_0  = input.c_0;      c_sct   = input.c_sct;    
 gam0 = input.gamma_0;  gam_sct = input.gamma_0*c_0/c_sct;  

% Transform Cartesian coordinates to spherical coordinates
  xS = input.xS;  rS = sqrt(xS(1)^2+xS(2)^2+xS(3)^2);
  X1 = input.X1;  X2 = input.X2;  X3 = input.X3; 
  R  = sqrt(X1.^2+X2.^2+X3.^2+1e-16);    % add small value to avoid zero R 

% (1) Compute incident wave in closed form --------------------------------
      DIS   = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
      u_inc = exp(-gam0*DIS)./(4*pi*DIS);
  
% (2) Compute coefficients of series expansion ----------------------------
      N = 50;                               % increase N for more accuracy 
      A = zeros(1,N+1);  B = zeros(1,N+1);
      arg0 = gam0 * input.a;  args = gam_sct * input.a;
      for n = 0 : N  
       Ib0 = besseli(n+1/2,arg0); dIb0 = besseli(n+3/2,arg0) + n/arg0*Ib0; 
       Ibs = besseli(n+1/2,args); dIbs = besseli(n+3/2,args) + n/args*Ibs; 
       Kb0 = besselk(n+1/2,arg0); dKb0 =-besselk(n+3/2,arg0) + n/arg0*Kb0;
       A(n+1) = - (gam_sct *dIbs*Ib0 - gam0 *dIb0*Ibs)  ...
                 /(gam_sct *dIbs*Kb0 - gam0 *dKb0*Ibs);
   % factor sqrt(pi^2/4/arg0/args) in numerator and denominator is omitted
       B(n+1) = sqrt(args/arg0)* (Ib0 + A(n+1)*Kb0) / Ibs;  
      end

% (3) Compute exterior, interior and total fields -------------------------
      DIS = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
      COS = (R.^2 + rS^2 - DIS.^2) ./ (2*rS*R);         %  cos(x,xS) 
      Pn = zeros(size(R)); Pn_1 = Pn; Pn_2 = Pn;  
      u_rfl = zeros(size(R));  
      u_int = zeros(size(R));
      for n = 0 : N
       [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);
       arg = gam0*R;   fctr=sqrt(pi/2./arg); Kb0=fctr.*besselk(n+1/2,arg);
       arg = gam0*rS;  fctr=sqrt(pi/2 /arg); KbS=fctr *besselk(n+1/2,arg);
       arg = gam_sct*R;fctr=sqrt(pi/2./arg); Ibs=fctr.*besseli(n+1/2,arg);                  
       u_rfl = u_rfl + A(n+1)*(2*n+1) * Kb0 .* KbS .* Pn;  
       u_int = u_int + B(n+1)*(2*n+1) * Ibs .* KbS .* Pn;
      end  
      u_rfl(ceil(input.N1/2),ceil(input.N2/2),ceil(input.N3/2)) = 0;
                                                % Avoid NaN number at R=0
      u_rfl = u_rfl * gam0 / (2*pi^2);  
      u_int = u_int * gam0 / (2*pi^2); 
  
      u = (u_rfl+u_inc) .* (R >= input.a) +  u_int .* (R < input.a);  
 
      plotWavefield(u_inc,u,input)
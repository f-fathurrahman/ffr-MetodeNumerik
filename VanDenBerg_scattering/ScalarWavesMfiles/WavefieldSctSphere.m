clear all; clc; close all;
input = init();
c_0  = input.c_0;       c_sct   = input.c_sct;       
gam0 = input.gamma_0;   gam_sct = input.gamma_0*c_0/c_sct;  

if exist(fullfile(cd, 'DATA3D.mat'), 'file');   delete DATA3D.mat;   end  

% (1) Compute coefficients of series expansion ----------------------------
      N = 50;                               % increase N for more accuracy 
      A = zeros(1,N+1);
      arg0 = gam0 * input.a;       args = gam_sct * input.a;
      for n = 0 : N  
       Ib0 = besseli(n+1/2,arg0); dIb0 = besseli(n+3/2,arg0) + n/arg0*Ib0; 
       Ibs = besseli(n+1/2,args); dIbs = besseli(n+3/2,args) + n/args*Ibs; 
       Kb0 = besselk(n+1/2,arg0); dKb0 =-besselk(n+3/2,arg0) + n/arg0*Kb0;
       A(n+1) =  -(gam_sct *dIbs*Ib0 - gam0 *dIb0*Ibs)  ... 
                / (gam_sct *dIbs*Kb0 - gam0 *dKb0*Ibs);
   % factor sqrt(pi^2/4/arg0/args) in numerator and denominator is omitted
      end

% (2) Compute reflected field at receivers (data) -------------------------
      xR  = input.xR;  rR = sqrt(xR(1,:).^2 + xR(2,:).^2 + xR(3,:).^2);
      xS  = input.xS; rS = sqrt(xS(1)^2+xS(2)^2+xS(3)^2); 
      DIS = sqrt((xR(1,:)-xS(1)).^2+(xR(2,:)-xS(2)).^2+(xR(3,:)-xS(3)).^2);
      COS = (rR.^2 + rS^2 - DIS.^2) ./ (2*rR*rS);         %  cos(xR,xS)
      Pn  = zeros(1,input.NR); Pn_1 = Pn; Pn_2 = Pn;
      data3D = zeros(1,input.NR); 
      for n = 0 : N
       [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);
       arg = gam0*rR; fctr=sqrt(pi/2./arg); Kb0=fctr.* besselk(n+1/2,arg);
       arg = gam0*rS; fctr=sqrt(pi/2 /arg); KbS=fctr * besselk(n+1/2,arg);
       data3D = data3D + A(n+1)*(2*n+1) * Kb0 .* KbS .* Pn;                     
      end % n_loop
      data3D = data3D * gam0 / (2*pi^2);  
                                                      
      displayData(data3D,input);                       save DATA3D data3D;
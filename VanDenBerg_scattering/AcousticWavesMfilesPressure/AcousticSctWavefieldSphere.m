clear all; clc; close all; clear workspace;    input = initAC();
c_0   = input.c_0;      c_sct   = input.c_sct;  
rho_0 = input.rho_0;    rho_sct = input.rho_sct; 
gam0  = input.gamma_0;  gam_sct = input.gamma_0 * c_0/c_sct; 
Z_0   = rho_0 * c_0;    Z_sct   = rho_sct * c_sct;

if exist(fullfile(cd, 'DATA3D.mat'), 'file');   delete DATA3D.mat;   end   
% (1) Compute coefficients of series expansion ----------------------------
      N = 50;                                 % increase N for more accuracy 
      A = zeros(N+1); B = zeros(N+1);
      arg0 = gam0    * input.a;  
      args = gam_sct * input.a;
    for n = 0 : N  
      Ib0 = besseli(n+1/2,arg0);  dIb0 = besseli(n+3/2,arg0) + n/arg0*Ib0; 
      Ibs = besseli(n+1/2,args);  dIbs = besseli(n+3/2,args) + n/args*Ibs; 
      Kb0 = besselk(n+1/2,arg0);  dKb0 =-besselk(n+3/2,arg0) + n/arg0*Kb0;
      denominator = (1/Z_sct) * dIbs*Kb0 - (1/Z_0) * dKb0*Ibs;
      A(n+1) = -((1/Z_sct) * dIbs*Ib0 - (1/Z_0) * dIb0*Ibs) / denominator; 
      B(n+1) = (Ib0 + A(n+1)*Kb0)/Ibs;
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
         factor        = (2*n+1)*besselk(n+1/2,gam0*rS).*Pn;
         data3D        = data3D + A(n+1).*factor.* besselk(n+1/2,gam0*rR);                     
      end % n_loop
      data3D = data3D ./ (4*pi*sqrt(rS*rR));  
  
      displayData(data3D,input);                       save DATA3D data3D;
clear all; clc; close all; clear workspace
input = initAC();                                        
a     = input.a;             xS = input.xS;
c_0   = input.c_0;      c_sct   = input.c_sct;  
rho_0 = input.rho_0;    rho_sct = input.rho_sct; 
gam0  = input.gamma_0;  gam_sct = input.gamma_0 * c_0/c_sct;
Z_0   = rho_0 * c_0;    Z_sct   = rho_sct * c_sct;

if exist(fullfile(cd, 'DATA1D.mat'), 'file');   delete DATA1D.mat;   end
% (1) Compute incident pressure wave field --------------------------------
      p_inc_a = 1/(2*gam0) * exp(gam0*(xS+a));
% (2) Compute coefficents A and B of internal field -----------------------
      Rho = (1/Z_0-1/Z_sct) / (1/Z_0+1/Z_sct); 
      Tau =  2/Z_0          / (1/Z_0+1/Z_sct);
      denominator = 1 - Rho * Rho * exp(-4*gam_sct*a); 
      A = Tau / denominator;  
      B = - Rho * A * exp(-2*gam_sct*a);
% (3) Compute reflection and transmision factors --------------------------
      R = A + exp(-2*gam_sct*a) * B - 1;  
      T = (B + exp(-2*gam_sct*a) * A) * exp(2*gam0*a);  
% (4) Compute receiver data for reflection and transmission ---------------
      xR = input.xR(1,1);    
      data1D(1) = p_inc_a * R * exp(gam0*(xR+a));
      xR = input.xR(1,2);
      data1D(2) = p_inc_a * (T-1) * exp(-gam0*(xR+a));    
      displayData(data1D,input);                       save DATA1D data1D;
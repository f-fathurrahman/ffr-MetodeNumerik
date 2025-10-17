clear all; clc; close all;
input = init();       xS = input.xS; a = input.a;                                       
gam0 = input.gamma_0; gam_sct = input.gamma_0 * input.c_0 / input.c_sct; 

if exist(fullfile(cd, 'DATA1D.mat'), 'file');   delete DATA1D.mat;   end
% (1) Compute incident wave -----------------------------------------------
      u_inc_a = 1/(2*gam0) * exp(gam0*(xS+a));
% (2) Compute coefficents A and B of internal field -----------------------
      rho         = (gam0-gam_sct) / (gam0+gam_sct); 
      tau         =  2*gam0        / (gam0+gam_sct);
      denominator = 1 - rho * rho * exp(-4*gam_sct*a); 
      A           = tau / denominator;  
      B           = - rho * A * exp(-2*gam_sct*a);
% (3) Compute reflection and transmision factors --------------------------
      R =  A + exp(-2*gam_sct*a) * B - 1;  
      T =  (B + exp(-2*gam_sct*a) * A) * exp(2*gam0*a);  
% (4) Compute receiver data for reflection and transmission ---------------
             xR = input.xR(1,1);    
      data1D(1) = u_inc_a * R * exp(gam0*(xR+a));
             xR = input.xR(1,2);
      data1D(2) = u_inc_a * (T-1) * exp(-gam0*(xR+a));   
      
      displayData(data1D,input);                       save DATA1D data1D;
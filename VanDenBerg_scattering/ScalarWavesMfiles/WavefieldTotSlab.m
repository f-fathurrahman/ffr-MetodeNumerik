clear all; clc; close all;
input = init();         X1 = input.X1; xS = input.xS;   a = input.a; 
gam0  = input.gamma_0;  gam_sct = input.gamma_0 * input.c_0 / input.c_sct; 

% (1) Compute incident wave -----------------------------------------------
      u_inc   = 1/(2*gam0) * exp(-gam0*abs(xS-X1));
      u_inc_a = 1/(2*gam0) * exp(gam0*(xS+a));
  
% (2) Compute coefficients A and B of internal field ----------------------
      rho         = (gam0-gam_sct) / (gam0+gam_sct); 
      tau         =  2*gam0        / (gam0+gam_sct);
      denominator = 1 - rho * rho * exp(-4*gam_sct*a); 
      A           = tau / denominator;  
      B           = - rho * A * exp(-2*gam_sct*a);
    
% (3) Compute reflection and transmision factors --------------------------
      R =  A + exp(-2*gam_sct*a) * B - 1;  
      T = (B + exp(-2*gam_sct*a) * A) * exp(2*gam0*a) ;

% (4) Compute total wave field distribution ------------------------------- 
      u = zeros(input.N1,1); 
      for i = 1:input.N1 
        if X1(i) < -a;  u(i) = u_inc(i)+R*u_inc_a*exp(gam0*(X1(i)+a)); end
        if abs(X1(i)) < a     
            u(i) = u_inc_a * ( A * exp(-gam_sct*(X1(i)+a)) ...
                                           + B * exp(gam_sct*(X1(i)-a)) );
        end
        if  X1(i) > a;  u(i) = T * u_inc_a * exp(-gam0.*(X1(i)+a)); end
      end % i-loop 
      plotWavefield(u_inc,u,input)
clear all; clc; close all;
input = init();  gam0 = input.gamma_0;    a = input.a;      
   
% Transform Cartesian coordinates to polar ccordinates  
  xS = input.xS; 
  rS = sqrt(xS(1)^2+xS(2)^2);  phiS  = atan2(xS(2),xS(1)); 
  X1 = input.X1; X2 = input.X2;  
  R  = sqrt(X1.^2+X2.^2);  PHI = atan2(X2,X1); 

% (1) Compute incident wave in closed form --------------------------------
      DIS         = sqrt(rS^2+R.^2-2*rS*R.*cos(phiS-PHI));
      u_inc_exact = 1/(2*pi) .* besselk(0,gam0*DIS);

% (2) Compute incident wave as Bessel series with -M:M terms --------------
      M = 50;                               % increase M for more accuracy
      u_inc = besselk(0,gam0*rS) .* besseli(0,gam0*R);   % zero order term
      for m = 1 : M           
        u_inc = u_inc + 2 * besselk(m,gam0*rS) ...
                              .* besseli(m,gam0*R) .* cos(m*(phiS-PHI));
      end % m_loop
      u_inc = 1/(2*pi) * u_inc;

% (3) Determine mean error and plot error in domain -----------------------   
      Error = u_inc - u_inc_exact;
      disp(['normalized error = '  ...
                       num2str(norm(Error(:),1)/norm(u_inc_exact(:),1))]);        
      set(figure,'Units','centimeters','Position',[5 5 18 10]);            
      subplot(1,3,1)
       IMAGESC(X1(:,1),X2(1,:),abs(u_inc_exact))
       title('\fontsize{10} 2D: abs(u^{inc}_{exact})');
       hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');
      subplot(1,3,2)
       IMAGESC(X1(:,1),X2(1,:),abs(u_inc)) 
       title('\fontsize{10} 2D: abs(u^{inc})');  
       hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');
     subplot(1,3,3)
      IMAGESC(X1(:,1),X2(1,:),abs(Error))
      title('\fontsize{10} 2D: abs(Error)');   
      hold on; phi = 0:.01:2*pi;  plot(a*cos(phi),a*sin(phi),'w');
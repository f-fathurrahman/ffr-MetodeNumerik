clear all; clc; close all;
input = init(); 
gam0 = input.gamma_0;   a = input.a;   
   
% Transform Cartesian coordinates to spherical coordinates
  xS = input.xS;   
  rS = sqrt(xS(1)^2+xS(2)^2+xS(3)^2);
  X1 = input.X1;   X2 = input.X2;   X3 = input.X3; 
  R  = sqrt(X1.^2+X2.^2+X3.^2+1e-16);    % add small value to avoid zero R  

% (1) Compute incident wave in closed form --------------------------------
      DIS         = sqrt( (X1-xS(1)).^2 + (X2-xS(2)).^2 + (X3-xS(3)).^2 );
      u_inc_exact = exp(-gam0*DIS)./(4*pi*DIS);
    
% (2) Compute incident wave as Bessel series with 0:N terms ---------------
      N = 50;                               % increase N for more accuracy
      COS   = (R.^2 + rS^2 - DIS.^2) ./ (2*rS*R);              % cos(x,xS)
      Pn    = zeros(size(R)); Pn_1 = Pn; Pn_2 = Pn;
      u_inc = zeros(size(R));  
      for n = 0 : N                           
        [Pn,Pn_1,Pn_2] = Legendre(n,COS,Pn,Pn_1,Pn_2);
        arg = gam0*R; fctr=sqrt(pi/2./arg); Ib0=fctr.* besseli(n+1/2,arg);
        arg = gam0*rS;fctr=sqrt(pi/2 /arg); KbS=fctr * besselk(n+1/2,arg);
        u_inc = u_inc + (2*n+1) * Ib0 .* KbS .* Pn;        
      end % n_loop
        u_inc = u_inc * gam0 / (2*pi^2);
  
% (3) Determine and plot relative error in domain ------------------------- 
      Error = u_inc - u_inc_exact;
      disp(['mean_square_error = '   ...
                       num2str(norm(Error(:),1)/norm(u_inc_exact(:),1))]);
    % plot at  cross-section at either x3 = 0 or x3 = dx/2 
      N1 = input.N1;  N2 = input.N2;  N3 = input.N3;
      N_cross =floor(N3/2+1); 
      set(figure,'Units','centimeters','Position',[5 5 18 10]);    
      subplot(1,3,1)
       matrix2D = reshape(abs(u_inc_exact(:,:, N_cross)),N1,N2); 
       IMAGESC(X1(:,1,1),X2(1,:,1),matrix2D)
       title('\fontsize{11} 3D: abs(u^{inc}_{exact})');
       hold on; phi = 0:.01:2*pi;  plot(a*cos(phi),a*sin(phi),'w');
      subplot(1,3,2) 
       matrix2D = reshape(abs(u_inc(:,:, N_cross)),N1,N2);
       IMAGESC(X1(:,1,1),X2(1,:,1),matrix2D); 
       title('\fontsize{11} 3D: abs(u^{inc})');
       hold on; phi = 0:.01:2*pi;  plot(a*cos(phi),a*sin(phi),'w'); 
      subplot(1,3,3)
       matrix2D = reshape(abs(Error(:,:, N_cross)),N1,N2);
       IMAGESC(X1(:,1,1),X2(1,:,1),matrix2D)
       title('\fontsize{11} 3D: abs({Error})');
       hold on; phi = 0:.01:2*pi;  plot(a*cos(phi),a*sin(phi),'w');
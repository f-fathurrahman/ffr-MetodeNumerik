function [input] = initAcousticContrast(input)

x1 = input.X1(:,1);  x2 = input.X2(1,:);   dx = input.dx;
N1 = input.N1;       N2 = input.N2;

input.a  = 40;  % half width slab / radius circle cylinder / radius sphere
input.xO(1) = input.a / 2;   % center coordinates of circle
input.xO(2) = input.a / 3; 
          R = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2);
shape = (R < input.a);

% (1) Compute compressibbily contrast and mass density contrast ----------
      kappa_0        = 1 / (input.rho_0   * input.c_0^2);
      kappa_sct      = 1 / (input.rho_sct * input.c_sct^2);
      contrast_kappa = 1 - kappa_sct /kappa_0;
      input.CHI_kap  = contrast_kappa .* shape;   
      contrast_rho   = 1 - input.rho_sct / input.rho_0; 
      input.CHI_rho  = contrast_rho .* shape;

   % Plot contrast values
     set(figure(1),'Units','centimeters','Position',[5 1 18 12]);  
     subplot(1,2,1);
       IMAGESC(x1,x2, real(input.CHI_kap));  
       title('\fontsize{13}\chi^{\kappa} = 1 - \kappa_{sct}/ \kappa_{0}');
     subplot(1,2,2);      
       IMAGESC(x1,x2,input.CHI_rho);  
       title('\fontsize{13}\chi^{\rho} = 1 - \rho_{sct}/ \rho_{0}');
  
% (2) Compute volume contrast and interface contrast ---------------------
      c_contrast = 1 - input.c_0^2/input.c_sct^2;  
      input.CHI  = c_contrast .* shape;    
      rho        = input.rho_sct .* shape + (1-shape) .* input.rho_0;

      Rfl{1} = zeros(N1,N2);   Rfl{2} = zeros(N1,N2);
      Rfl{1}(1:N1-1,:) = (rho(2:N1,:) - rho(1:N1-1,:)) ...
                                      ./(rho(2:N1,:) + rho(1:N1-1,:));   
      Rfl{2}(:,1:N2-1) = (rho(:,2:N2) - rho(:,1:N2-1)) ...
                                      ./(rho(:,2:N2) + rho(:,1:N2-1));                        
      input.Rfl = Rfl;                      

   % Plot contrast values
     set(figure(2),'Units','centimeters','Position',[5 1 18 12]); 
     subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2,Rfl{1});
       axis('ij','equal','tight');  title('\fontsize{13} R^{(1)}');
       colorbar('hor')
    subplot(1,2,2);  
      IMAGESC(x1,x2+dx/2,Rfl{2});
      axis('ij','equal','tight');  title('\fontsize{13} R^{(2)}');
      colorbar('hor');  colormap(jet);  
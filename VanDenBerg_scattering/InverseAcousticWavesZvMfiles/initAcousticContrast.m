function [input] = initAcousticContrast(input)

input.a  = 40;  % half width slab / radius circle cylinder / radius sphere
   
input.xO(1) = input.a / 2;   % center coordinates of circle
input.xO(2) = input.a / 3; 
          R = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2);

shape = (R < input.a);

% (1) Compute compressibbily contrast (kappa = 1/(rho*c^2) ---------------
      kappa_0        = 1 / (input.rho_0   * input.c_0^2);
      kappa_sct      = 1 / (input.rho_sct * input.c_sct^2);
      contrast_kappa = 1 - kappa_sct /kappa_0;
      input.CHI_kap  = contrast_kappa .* shape;   

% (2) Compute mass density contrast --------------------------------------
      contrast_rho   = 1 - input.rho_sct / input.rho_0; 
      input.CHI_rho  = contrast_rho .* shape;

  % Plot contrast values
    x1 = input.X1(:,1);   x2 = input.X2(1,:); 
    figure(1);
      IMAGESC(x1,x2, real(input.CHI_kap));  
      title('\fontsize{13} \chi^{\kappa} = 1 - \kappa_{sct}/ \kappa_{0}');
    figure(2); 
      mesh(abs(input.CHI_kap),'edgecolor', 'k'); 
      view(37.5,45); axis xy; axis('off'); axis('tight') 
      axis([1 input.N1 1 input.N2  0 10000*max(abs(input.CHI_kap(:)))])    
    figure(3);  
      IMAGESC(x1,x2,input.CHI_rho);  
      title('\fontsize{13}\chi^{\rho} = 1 - \rho_{sct}/ \rho_{0}');
    figure(4); 
      mesh(abs(input.CHI_rho),'edgecolor', 'k'); 
      view(37.5,45); axis xy; axis('off'); axis('tight') 
      axis([1 input.N1 1 input.N2  0 1.5*max(abs(input.CHI_rho(:)))])
function [CHI,r_D,eta_D,ErrorD] = ...
                          UpdateContrastGaussNewton(w,u,u_inc,input)
global it; NS = input.NS; 

% (1) Compute CHI_csi and determine ErrorD
      INuw = zeros(input.N1,input.N2);    INuu = zeros(input.N1,input.N2);  
      for q = 1 : input.NS     
        INuw = INuw + (conj(u{q}) .* w{q}); 
        INuu = INuu + real(conj(u{q}) .* u{q});
      end %q_loop  
      CHI_csi = real(INuw./INuu); 

% (2) Gauss Newton method for unknown location, contrast and radius
      [CHI,input] = GaussNewtonCSI(CHI_csi,input);
      ErrorCHI = norm(CHI(:)-input.CHI(:));
      txt = ['         ErrorCHI = ', num2str(ErrorCHI/norm(CHI(:)))];
      disp(['Iteration = ', num2str(it),txt]);
    
      Norm_D = 0;   Norm = 0;    r_D{q} = cell(1,input.NS); 
      for q = 1 : NS
        r_D{q} = CHI .* u{q} - w{q};  
        Norm_D = Norm_D + norm(r_D{q}(:))^2;
        Norm   = Norm + norm(CHI(:).*u_inc{q}(:))^2;
      end %q_loop 
      eta_D = 1 / (Norm); 
      ErrorD = eta_D * Norm_D;  

% (3) Show intermediate pictures of reconstruction  
      x1 = input.X1(:,1);   x2 = input.X2(1,:);
      figure(5); IMAGESC(x1,x2,real(CHI)); 
      title('\fontsize{13} Reconstructed Contrast Re[\chi^{anl}] ');
      pause(0.1)  
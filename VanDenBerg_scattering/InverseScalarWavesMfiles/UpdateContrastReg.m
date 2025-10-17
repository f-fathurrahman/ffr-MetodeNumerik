function [CHI,r_D,eta_D,ErrorD] = ...
                          UpdateContrastReg(w,u,u_inc,ErrorM,ErrorD,input)

% (1) Compute CHI_csi 
      INuw = zeros(input.N1,input.N2);   INuu = zeros(input.N1,input.N2);  
      for q = 1 : input.NS     
        INuw = INuw + conj(u{q}) .* w{q}; 
        INuu = INuu + conj(u{q}) .* u{q};
      end %q_loop  
      CHI = INuw./INuu; 

% (2) Add regularization by Jacobi iteration of Euler-Lagrange equation
      CHI = TVregularizer(CHI,ErrorM,ErrorD,input);
        
      % Enhance TV minimization by repeating the statements:
          % CHI = TVregularizer(CHI,ErrorM,ErrorD,input);  
          % CHI = TVregularizer(CHI,ErrorM,ErrorD,input);

% (3) Update the object error in D
      Norm_D = 0;   Norm = 0;   
      r_D = cell(1,input.NS); 
      for q = 1 : input.NS
        r_D{q} = CHI .* u{q} - w{q};  
         Norm_D = Norm_D + norm(r_D{q}(:))^2;
         Norm   = Norm + norm(CHI(:).*u_inc{q}(:))^2;
      end %q_loop 
      eta_D = 1 / (Norm); 
      ErrorD = eta_D * Norm_D;  
      disp([' ErrorM = ',num2str(ErrorM),'   ErrorD = ',num2str(ErrorD)]);
      
% (4) Show intermediate pictures of reconstruction 
      x1 = input.X1(:,1);   
      x2 = input.X2(1,:);
      figure(5); 
      IMAGESC(x1,x2,real(CHI)); 
      title('\fontsize{13} Reconstructed Contrast Re[\chi^c] ');
      pause(0.1) 
function [input] = initAcousticContrastIntf(input)
global nDIM; 

if     nDIM == 1;  R = sqrt(input.X1.^2);               [N1,~]  = size(R);   
elseif nDIM == 2;  R = sqrt(input.X1.^2 + input.X2.^2); [N1,N2] = size(R);
elseif nDIM == 3;  R = sqrt(input.X1.^2 + input.X2.^2 + input.X3.^2);
                   [N1,N2,N3] = size(R); 
end % if

% (1) Compute wave speed contrast mass density contrast -------------------
               a = input.a;
      c_contrast = 1 - input.c_0^2/input.c_sct^2;  
      input.CHI  = c_contrast * (R < a);   

% (2) Compute mass density and 'reflection factors' -----------------------
      rho = input.rho_sct * (R < a) + (R >= a) * input.rho_0; 

      if nDIM == 1 
              Rfl{1} = zeros(N1,1); 
      Rfl{1}(1:N1-1) = (rho(2:N1)-rho(1:N1-1))./(rho(2:N1)+rho(1:N1-1)); 
      elseif nDIM == 2 
                  Rfl{1} = zeros(N1,N2); 
                  Rfl{2} = zeros(N1,N2);
        Rfl{1}(1:N1-1,:) = (rho(2:N1,:) - rho(1:N1-1,:)) ...
                         ./(rho(2:N1,:) + rho(1:N1-1,:));   
        Rfl{2}(:,1:N2-1) = (rho(:,2:N2) - rho(:,1:N2-1)) ...
                         ./(rho(:,2:N2) + rho(:,1:N2-1));     
      elseif nDIM == 3 
                  Rfl{1} = zeros(N1,N2,N3); 
                  Rfl{2} = zeros(N1,N2,N3);
                  Rfl{3} = zeros(N1,N2,N3); 
      Rfl{1}(1:N1-1,:,:) = (rho(2:N1,:,:) - rho(1:N1-1,:,:)) ...
                           ./(rho(2:N1,:,:) + rho(1:N1-1,:,:));   
      Rfl{2}(:,1:N2-1,:) = (rho(:,2:N2,:) - rho(:,1:N2-1,:)) ...
                           ./(rho(:,2:N2,:) + rho(:,1:N2-1,:));
      Rfl{3}(:,:,1:N3-1) = (rho(:,:,2:N3) - rho(:,:,1:N3-1)) ...
                           ./(rho(:,:,2:N3) + rho(:,:,1:N3-1));                      
     end % if

     input.Rfl = Rfl;                   clear input.CHI_rho  input.CHI_kap
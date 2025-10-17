function [input] = initEMcontrastIntf(input)
global nDIM; 
  
if     nDIM == 2;  R = sqrt(input.X1.^2 + input.X2.^2); [N1,N2] = size(R);
elseif nDIM == 3;  R = sqrt(input.X1.^2 + input.X2.^2 + input.X3.^2);
                   [N1,N2,N3] = size(R); 
end % if

% Compute permittivity type of 'reflection factors' ----------------------
   a = input.a;
  mu = input.mu_sct * (R < a) + (R >= a) * 1; 


if nDIM == 2; 
              Rfl{1} = zeros(N1,N2); 
              Rfl{2} = zeros(N1,N2);
    Rfl{1}(1:N1-1,:) = (mu(2:N1,:) - mu(1:N1-1,:)) ...
                         ./(mu(2:N1,:) + mu(1:N1-1,:));   
    Rfl{2}(:,1:N2-1) = (mu(:,2:N2) - mu(:,1:N2-1)) ...
                         ./(mu(:,2:N2) + mu(:,1:N2-1));     
elseif nDIM == 3; 
              Rfl{1} = zeros(N1,N2,N3); 
              Rfl{2} = zeros(N1,N2,N3);
              Rfl{3} = zeros(N1,N2,N3); 
  Rfl{1}(1:N1-1,:,:) = (mu(2:N1,:,:) - mu(1:N1-1,:,:)) ...
                           ./(mu(2:N1,:,:) + mu(1:N1-1,:,:));   
  Rfl{2}(:,1:N2-1,:) = (mu(:,2:N2,:) - mu(:,1:N2-1,:)) ...
                           ./(mu(:,2:N2,:) + mu(:,1:N2-1,:));
  Rfl{3}(:,:,1:N3-1) = (mu(:,:,2:N3) - mu(:,:,1:N3-1)) ...
                           ./(mu(:,:,2:N3) + mu(:,:,1:N3-1));                      
end % if

input.Rfl = Rfl;                                           
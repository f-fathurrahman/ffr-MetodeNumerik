function [input] = initEMcontrastIntf(input) 
 
R = sqrt(input.X1.^2 + input.X2.^2); 
[N1,N2] = size(R);
       
% Compute permittivity type of 'reflection factors' ----------------------
    a = input.a;
  eps = input.eps_sct * (R < a) + (R >= a) * 1; 
    
  Rfl{1} = zeros(N1,N2); 
  Rfl{2} = zeros(N1,N2);
  Rfl{1}(1:N1-1,:) = (eps(2:N1,:) - eps(1:N1-1,:)) ...
                                 ./(eps(2:N1,:) + eps(1:N1-1,:));   
  Rfl{2}(:,1:N2-1) = (eps(:,2:N2) - eps(:,1:N2-1)) ...
                                 ./(eps(:,2:N2) + eps(:,1:N2-1));     
 input.Rfl = Rfl;                                           
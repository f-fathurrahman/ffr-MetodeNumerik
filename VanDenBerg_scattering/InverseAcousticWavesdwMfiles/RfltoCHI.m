function [CHI] = RfltoCHI(Rfl,input)

N1 = input.N1;  N2 = input.N2;     
  
% Determine the mass-density distribution 
  factor1 = (1 + Rfl{1}) ./ (1 - Rfl{1});
% and integrate in  forward x_1 direction
  rhoF1 = ones(N1,N2);  
  for i = 1 : N1-1
    rhoF1(i+1,:) = rhoF1(i,:) .* factor1(i,:);
  end 
% and integrate in backward x_1 direction
  rhoB1 = ones(N1,N2);  
  for i = 2 : N1
      I = N1-i+1;    rhoB1(I,:) = rhoB1(I+1,:) ./ factor1(I+1,:);
  end 
  i = 2:N1-1; rhoB1(i,:) = rhoB1(i-1,:);          % with shift correction
  
% Determine the mass-density distribution 
  factor2 = (1 + Rfl{2}) ./ (1 - Rfl{2});  
% and integrate in forward x_2 direction  
  rhoF2 = ones(N1,N2);  
  for j = 1 : N2-1
    rhoF2(:,j+1) = rhoF2(:,j) .* factor2(:,j);
  end
% and integrate in backward x_2 direction
  rhoB2 = ones(N1,N2);  
  for j = 2 : N2
      J = N2-j+1;    rhoB2(:,J) = rhoB2(:,J+1) ./ factor2(:,J+1);
  end   
  j = 2:N2-1; rhoB2(:,j) = rhoB2(:,j-1);          % with shift correction 

% determine mean value and mass-density contrast CHI
  CHI = 1- (rhoF1 + rhoB1 + rhoF2 + rhoB2)/4;
 
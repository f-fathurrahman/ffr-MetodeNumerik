function [CHI] = TVregularizer(CHI,ErrorM,ErrorD,input)
N1 = input.N1;  N2 = input.N2;

% (1) Determine weights b^2 
  Ext_CHI                = zeros(N1+2,N2+2);   
  Ext_CHI(2:N1+1,2:N2+1) = CHI(1:N1,1:N2);
  
  Fwd = zeros(N1+2,N2+2);  
  Bwd = zeros(N1+2,N2+2);
  
  % Symmetric forward and backward finite differences
  i = 2:N1+1;    %-----------------------------------------------
     j = 1:N2+1;  Fwd(i,j) = abs(Ext_CHI(i,j+1)-Ext_CHI(i,j)).^2;
     j = 2:N2+2;  Bwd(i,j) = abs(Ext_CHI(i,j)-Ext_CHI(i,j-1)).^2;
  j = 2:N2+1;    %--------------------------------------------------------
     i = 1:N1+1;  Fwd(i,j) = Fwd(i,j)+abs(Ext_CHI(i+1,j)-Ext_CHI(i,j)).^2;  
     i = 2:N1+2;  Bwd(i,j) = Bwd(i,j)+abs(Ext_CHI(i,j)-Ext_CHI(i-1,j)).^2;

  % Take delta^2 either as ErrorD or as mean value of profile gradient
     deltan = ErrorD;
    %deltan = mean(Fwd(:)+Bwd(:))/2;

  i = 2:N1+1; j = 2:N2+1; %-------------------------------
      bn(i-1,j-1) = 1 ./ ((Bwd(i,j)+Fwd(i,j))/2 + deltan);
    
% (2) Determine regularized CHI by forward and backward differences
  Ext_bn                = zeros(N1+2,N2+2);    
  Ext_bn(2:N1+1,2:N2+1) = bn(1:N1,1:N2); 
  Ext_CHI               = zeros(N1+2,N2+2);   
  Ext_CHI(2:N1+1,2:N2+1)= CHI(1:N1,1:N2);
  
  D     = zeros(N1+2,N2+2);  
  R_CHI = zeros(N1+2,N2+2);  
  
  factor  = max(ErrorM,0.01);  % strength of regularization
  MeanCHI = mean(abs(CHI(:)).^2);
  
% (3) UPDATE THE CONTRAST using a single Jacobi iteration 
  i = 2:N1+1; j = 2:N2+1; %-----------------------------------------------
    D(i,j) = 1 + factor * MeanCHI /2                                   ...
                   .*( Ext_bn(i+1,j)+2*Ext_bn(i,j)+Ext_bn(i,j+1)       ... 
                     + Ext_bn(i-1,j)+2*Ext_bn(i,j)+Ext_bn(i,j-1) );
                 
    R_CHI(i,j) =  - factor * MeanCHI /2                                ...
                   .*( (Ext_bn(i+1,j)+Ext_bn(i,j)) .* Ext_CHI(i+1,j)   ...
                      +(Ext_bn(i-1,j)+Ext_bn(i,j)) .* Ext_CHI(i-1,j)   ...    
                      +(Ext_bn(i,j+1)+Ext_bn(i,j)) .* Ext_CHI(i,j+1)   ...
                      +(Ext_bn(i,j-1)+Ext_bn(i,j)) .* Ext_CHI(i,j-1) );              
    CHI(i-1,j-1) = (1./D(i,j)).*(Ext_CHI(i,j)-R_CHI(i,j));
function [CHI_eps,r_D,eta_D,ErrorD] = ...
                          UpdateEMContrast(wE,E,Einc,ErrorM,ErrorD,input) 

% Determine permittivity contrast
  INE_wE = zeros(input.N1,input.N2);     INE_E = zeros(input.N1,input.N2);  
  for q = 1 : input.NS     
    INE_wE = INE_wE + conj(E{1,q}) .* wE{1,q} + conj(E{2,q}) .* wE{2,q}; 
    INE_E  = INE_E  + conj(E{1,q}) .*  E{1,q} + conj(E{2,q}) .*  E{2,q};
  end %q_loop  
  CHI_eps = (INE_wE) ./ real(INE_E); 

  if ErrorD < 1 
     CHI_eps = TVregularizer(CHI_eps,ErrorM,ErrorD,input); 
  end
 
% Compute residual error in object domain   
  Norm_D = 0;   Norm = 0;    r_D{q} = cell(2,input.NS); 
  for q = 1 : input.NS
    r_D{1,q} = CHI_eps .* E{1,q} - wE{1,q}; 
    r_D{2,q} = CHI_eps .* E{2,q} - wE{2,q};
    Norm_D   = Norm_D + norm(r_D{1,q}(:))^2 + norm(r_D{2,q}(:))^2 ;
    Norm     = Norm   + norm(CHI_eps(:) .* Einc{1,q}(:))^2  ...
                      + norm(CHI_eps(:) .* Einc{2,q}(:))^2;
  end %q_loop 
  eta_D = 1 / (Norm);  
  ErrorD = eta_D * Norm_D;  
  disp(['  ErrorM = ',num2str(ErrorM), '    ErrorD  = ',num2str(ErrorD)]);

% (3) Show intermediate pictures of reconstruction 
  x1 = input.X1(:,1);   x2 = input.X2(1,:);
  figure(5); IMAGESC(x1,x2,real(CHI_eps)); 
  title('\fontsize{13} Reconstructed Contrast Re[\chi^\epsilon] ');
  pause(0.1) 
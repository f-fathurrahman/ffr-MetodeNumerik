function [CHI_rho,rZv,eta_D,ErrorD] ...
                    = UpdateContrastwZv(wZv,Zv,Zvinc,ErrorM,ErrorD,input)
global nDIM  

rZv  = cell(nDIM,input.NS);          
% Determine mass-density contrast  
  INZv_wZv = zeros(input.N1,input.N2);     
  INZv_Zv  = zeros(input.N1,input.N2);  
  for q = 1 : input.NS     
    INZv_wZv = INZv_wZv + conj(Zv{1,q}) .* wZv{1,q}   ...
                        + conj(Zv{2,q}) .* wZv{2,q};
    INZv_Zv  = INZv_Zv  + conj(Zv{1,q}) .*  Zv{1,q}   ...
                        + conj(Zv{2,q}) .*  Zv{2,q};
  end %q_loop 
  CHI_rho  = INZv_wZv ./ INZv_Zv;    

  if ErrorD < 1 
     CHI_rho  = TVregularizer(CHI_rho,ErrorM,ErrorD,input); 
  end

% Compute residual error in object domain 
  Norm_D = 0;   Norm = 0;
  for q = 1 : input.NS
    rZv{1,q} = CHI_rho .* Zv{1,q} - wZv{1,q};
    rZv{2,q} = CHI_rho .* Zv{2,q} - wZv{2,q};
    Norm_D   = Norm_D + norm(rZv{1,q}(:))^2 + norm(rZv{2,q}(:))^2;
    Norm     = Norm   + norm(CHI_rho(:) .* Zvinc{1,q}(:))^2 ...
                      + norm(CHI_rho(:) .* Zvinc{2,q}(:))^2; 
  end %q_loop 
  eta_D = 1 / Norm;
  ErrorD = eta_D * Norm_D;
  disp(['  ErrorM = ',num2str(ErrorM), '    ErrorD  = ',num2str(ErrorD)]);

% Show intermediate pictures of reconstruction of compressibility contrast 
  x1 = input.X1(:,1);   x2 = input.X2(1,:); 
  figure(5); 
     IMAGESC(x1,x2,real(CHI_rho));  
     title('\fontsize{13} Reconstructed Contrast \chi^{\rho}'); 
  figure(6); 
     mesh(real(CHI_rho),'edgecolor', 'k'); 
     view(37.5,45);  axis('off');  
     axis([1 input.N1 1 input.N2 [-1 1.01]*max(real(CHI_rho(:)))])     
pause(0.1);    
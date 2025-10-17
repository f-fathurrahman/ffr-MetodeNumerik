function [wZv,vZv,ZvD,gZv,AN,r_M] ...
      = InvITERwZv(CHI_rho,dG_R,wZv,vZv,ZvD,g_1Zv,AN_1,r_M,r_D,  ...
                                                      eta_M,eta_D,input)
global nDIM it
LZv = cell(nDIM,input.NS); 

% determine conjugate gradient directions
  dummy=cell(1,nDIM); 
  for n=1:nDIM; dummy{n}=conj(CHI_rho).*r_D{n}; end
  [KZv]  = AdjKOPwZv(dummy,input); 
  [gZv] = AdjDOPMwZv(dG_R,r_M,input);  
      AN = 0;
  for n = 1:nDIM 
     gZv{n} = eta_M * gZv{n} + eta_D * (r_D{n} - KZv{n});    
         AN = AN + norm(gZv{n}(:))^2;
  end
  if it == 1   
      for n = 1 : nDIM;  vZv{n} = gZv{n};  end                  
  else
      BN = 0;
      for n = 1:nDIM;   BN = BN + sum(gZv{n}(:).*conj(g_1Zv{n}(:))); end
      gamma = real((AN-BN)/AN_1);  
      for n = 1 : nDIM; vZv{n} = gZv{n} + gamma * vZv{n};
      end
  end    
  
% determine step length alpha
   GMp   = DOPMwZv(dG_R,vZv,input);
   BN_M  = norm(GMp(:))^2;        
   [KZv] = KOPwZv(vZv,input);                    
   AN = 0; BN_D = 0;  
   for n = 1:nDIM   
     LZv{n} = vZv{n} - CHI_rho .* KZv{n};  
       AN   = AN + sum(gZv{n}(:).*conj(vZv{n}(:))); 
       BN_D = BN_D + norm(LZv{n}(:))^2;                                   
   end
   alpha = real(AN) / (eta_M*BN_M + eta_D*BN_D); 
   
 % update contrast sources    
   for n = 1:nDIM  
        wZv{n} = wZv{n} + alpha * vZv{n};
        ZvD{n} = ZvD{n} + alpha * KZv{n};
   end
%  update residual error in M
   r_M  = r_M - alpha * GMp(:);
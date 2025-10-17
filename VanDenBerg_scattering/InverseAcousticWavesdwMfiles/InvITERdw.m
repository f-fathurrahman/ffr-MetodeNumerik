 function [dw,vdw,P,g_dw,AN,r_M] ...
            = InvITERdw(Rfl,dw,vdw,P,g_1dw,AN_1,r_M,r_D,eta_M,eta_D,input)
 global nDIM it
 Ldv = cell(nDIM,input.NS); shape = cell(nDIM);
% determine conjugate gradient directions
  dummy = cell(1,nDIM);
  for n = 1:nDIM; dummy{n} = conj(Rfl{n}).*r_D{n}; end
  if  input.Kirchhoff == 0;  [Kdf] = AdjKOPdw(dummy,input); end
  if (input.Kirchhoff == 1)  || (input.Kirchhoff == 2)
                    Kdf{1}= 0*dummy{1}; Kdf{2}= 0*dummy{2}; end             
  [g_dw] = AdjDOPMdw(r_M,input);
      AN = 0;   
  for n = 1:nDIM 
     % set data gradients for small interface contrast to zero
     shape{n} = abs(Rfl{n}) > 0.10 * max(abs(Rfl{n}(:)));
     g_dw{n}  = shape{n} .* g_dw{n};
     g_dw{n}  = eta_M * g_dw{n} + eta_D * (r_D{n} - Kdf{n}); 
     AN = AN + norm(g_dw{n}(:))^2;
  end
  if it == 1  
      for n = 1:nDIM; vdw{n} = g_dw{n};  end                
  else 
      BN = 0;
      for n = 1:nDIM; BN = BN + sum(g_dw{n}(:).* conj(g_1dw{n}(:))); end
      gamma = real((AN-BN)/AN_1);   
      for n = 1:nDIM  
        % set previous directions for small interface contrast to zero 
        vdw{n} = g_dw{n} + gamma * shape{n} .* vdw{n}; 
      end
 end   
% determine step length alpha
   GMp   = DOPMdw(vdw,input);
   BN_M  = norm(GMp(:))^2; 
   AN_D = 0;
   BN_D = 0;
   if (input.Kirchhoff == 0)  || (input.Kirchhoff == 1)
      [Kdv] = KOPdw(vdw,input); end
   if input.Kirchhoff == 2;  Kdv{1}= 0*vdw{1};  Kdv{2}= 0*vdw{2}; end
   for n = 1:nDIM
     AN_D    = AN_D + sum(g_dw{n}(:).*conj(vdw{n}(:)));
    Ldv{n} = vdw{n} - Rfl{n} .* Kdv{n};  
     BN_D  = BN_D + norm(Ldv{n}(:))^2;  
   end
  alpha  = real(AN_D) / (eta_M*BN_M + eta_D*BN_D); 
 % update contrast sources and AN 
   for n = 1:nDIM
     dw{n} = dw{n} + alpha * vdw{n}; 
     P{n}  = P{n} + alpha * Kdv{n};
   end
 % update residual error in M
     r_M = r_M - alpha * GMp(:);
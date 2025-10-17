function [w_Zv] = ITERCGwZv(Zv_inc,input)
global nDIM;  CHI_rho = input.CHI_rho;

w_Zv=cell(1,nDIM); g_Zv=cell(1,nDIM); r_Zv=cell(1,nDIM); v_Zv=cell(1,nDIM);
itmax  = 200; Errcri = input.Errcri;  it = 0;  % initialization iteration
Norm_D = 0;
for n = 1:nDIM 
   w_Zv{n} = zeros(size(Zv_inc{n}));  
   r_Zv{n} = CHI_rho.*Zv_inc{n};     Norm_D = Norm_D + norm(r_Zv{n}(:))^2; 
end  
eta_D = 1 / Norm_D;                                % normalization factor 
Error = 1;   fprintf('Error =         %g',Error);

while (it < itmax) && ( Error > Errcri) 
% determine conjugate gradient direction v
  dummy=cell(1,nDIM); for n=1:nDIM; dummy{n}=conj(CHI_rho).*r_Zv{n}; end
  [KZv] = AdjKOPwZv(dummy,input); 
  AN = 0;
  for n = 1:nDIM 
   g_Zv{n} = r_Zv{n} - KZv{n};              % and window for small CHI_rho
   g_Zv{n} = (abs(CHI_rho) >= Errcri).*g_Zv{n}; 
   AN = AN + norm(g_Zv{n}(:))^2;
  end
  if it == 0  
      for n=1:nDIM; v_Zv{n} = g_Zv{n};  end                  
  else
      for n=1:nDIM; v_Zv{n} = g_Zv{n} + AN/AN_1 * v_Zv{n}; end
  end    
% determine step length alpha        
  [KZv] = KOPwZv(v_Zv,input); 
   BN = 0;   
   for n = 1:nDIM   
     KZv{n} = v_Zv{n} - CHI_rho .* KZv{n};   BN = BN + norm(KZv{n}(:))^2;                                   
   end
   alpha = AN / BN;
 % update contrast sources and AN   
   for n = 1:nDIM;  w_Zv{n} = w_Zv{n} + alpha * v_Zv{n}; end
   AN_1  = AN;
 % update residual errors 
   Norm_D = 0;
   for n = 1:nDIM   
     r_Zv{n} = r_Zv{n} - alpha * KZv{n};  
     Norm_D  = Norm_D + norm(r_Zv{n}(:))^2; 
   end    
   Error = sqrt(eta_D * Norm_D); fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
 end % CG iterations
 
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached:  err/norm = ' num2str(Error)]); end
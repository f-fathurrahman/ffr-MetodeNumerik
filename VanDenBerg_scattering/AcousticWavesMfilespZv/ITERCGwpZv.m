function[w_p,w_Zv] = ITERCGwpZv(p_inc,Zv_inc,input)
global nDIM;   CHI_kap = input.CHI_kap;   CHI_rho = input.CHI_rho;

w_Zv=cell(1,nDIM); g_Zv=cell(1,nDIM); r_Zv=cell(1,nDIM); v_Zv=cell(1,nDIM);

itmax = 1000;  Errcri = input.Errcri;  it = 0;  % initialization iteration
  w_p = zeros(size(p_inc)); r_p = CHI_kap.*p_inc; Norm0 = norm(r_p(:))^2;
for n = 1:nDIM
   w_Zv{n} = zeros(size(Zv_inc{n}));  
   r_Zv{n} = CHI_rho.*Zv_inc{n};      Norm0 = Norm0 + norm(r_Zv{n}(:))^2;
end 

CheckAdjointpZv(r_p,r_Zv,input);  % Check once ---------------------------
Error = 1;   fprintf('Error =         %g',Error);
while (it < itmax) && ( Error > Errcri) && (Norm0 > eps)
% determine conjugate gradient direction v
  dummy=cell(1,nDIM); for n = 1:nDIM; dummy{n}=conj(CHI_rho).*r_Zv{n}; end
  [Kp,KZv] = AdjKOPpZv(conj(CHI_kap).*r_p, dummy,input); 
       g_p =  r_p - Kp;                    %  and window for small CHI_kap
       g_p = (abs(CHI_kap) >= Errcri).*g_p;     AN = norm(g_p(:))^2;
  for n = 1:nDIM
   g_Zv{n} = r_Zv{n} - KZv{n};              % and window for small CHI_rho
   g_Zv{n} = (abs(CHI_rho) >= Errcri).*g_Zv{n}; AN = AN+norm(g_Zv{n}(:))^2;
  end
  if it == 0;  v_p=g_p; for n=1:nDIM; v_Zv{n}=g_Zv{n};  end                 
  else
   v_p=g_p+AN/AN_1*v_p; for n=1:nDIM; v_Zv{n}=g_Zv{n}+AN/AN_1*v_Zv{n}; end
  end   
% determine step length alpha        
  [Kp,KZv] = KOPpZv(v_p,v_Zv,input); 
   Kp = v_p - CHI_kap .* Kp;                 BN = norm(Kp(:))^2 ;   
   for n = 1:nDIM  
     KZv{n} = v_Zv{n} - CHI_rho .* KZv{n};   BN = BN + norm(KZv{n}(:))^2;                                   
   end
   alpha = AN / BN;
 % update contrast sources and AN   
   w_p   = w_p + alpha * v_p;   
   for n = 1:nDIM;  w_Zv{n} = w_Zv{n} + alpha * v_Zv{n}; end
   AN_1  = AN;
 % update residual errors 
   r_p  = r_p  - alpha * Kp;              Norm = norm(r_p(:))^2;
   for n = 1:nDIM   
     r_Zv{n} = r_Zv{n} - alpha * KZv{n};  Norm = Norm + norm(r_Zv{n}(:))^2; 
   end    
   Error = sqrt(Norm / Norm0); fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
end  % CG iterations
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached:  err/norm = ' num2str(Error)]); end
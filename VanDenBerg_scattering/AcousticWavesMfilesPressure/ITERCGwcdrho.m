function [w_c,w_drho]=ITERCGwcdrho(p_inc,input)
global nDIM; CHI_c = input.CHI_c;   CHI_drho = input.CHI_drho;
w_drho = cell(1,nDIM);   Kwdrho = cell(1,nDIM);  
r_drho = cell(1,nDIM);   g_drho = cell(1,nDIM);  v_drho = cell(1,nDIM);

itmax = 1000;  Errcri = input.Errcri;  it = 0;  % initialization iteration 
w_c = zeros(size(p_inc));  r_c = CHI_c.*p_inc;    Norm0=norm(r_c(:))^2;
for n = 1:nDIM 
  w_drho{n} = zeros(size(p_inc));  r_drho{n} = CHI_drho{n}.*p_inc; 
      Norm0 = Norm0 + norm(r_drho{n}(:))^2;
end
CheckAdjointwcdrho(r_c,r_drho,input);  % Check once-----------------------
Error = 1;   fprintf('Error =         %g',Error);   
while (it < itmax) && ( Error > Errcri) && (Norm0 > eps)
  dummy = cell(1,nDIM);        % determine conjugate gradient direction v
  for n = 1:nDIM; dummy{n} = conj(CHI_drho{n}).*r_drho{n}; end
   [Kc,Kdrho] = AdjKOPwcdrho(conj(CHI_c).*r_c, dummy,input); 
          g_c =  r_c - Kc;                  %  and window for small CHI_c
          g_c = (abs(CHI_c) >= Errcri).*g_c;    AN = norm(g_c(:))^2;
  for n = 1:nDIM 
   g_drho{n} = r_drho{n} - Kdrho{n};        % and window for small CHI_rho
   g_drho{n} = (abs(CHI_drho{n}) >= Errcri).*g_drho{n}; 
   AN = AN+norm(g_drho{n}(:))^2;
  end
  if it == 0; v_c=g_c; for n=1:nDIM; v_drho{n}=g_drho{n};  end                  
  else
  v_c=g_c+AN/AN_1*v_c; for n=1:nDIM; v_drho{n}=g_drho{n}+AN/AN_1*v_drho{n};
  end
  end    
% determine step length alpha        
  [Kw] = KOPwcdrho(v_c,v_drho,input); 
   Kwc = v_c - CHI_c .* Kw;                BN = norm(Kwc(:))^2;   
  for n = 1:nDIM   
   Kwdrho{n} = v_drho{n}-CHI_drho{n}.*Kw;  BN = BN + norm(Kwdrho{n}(:))^2;                                   
  end
   alpha = AN / BN;     
   w_c   = w_c + alpha * v_c;    % update contrast sources 
   for n = 1:nDIM;  w_drho{n} = w_drho{n} + alpha * v_drho{n}; end
   AN_1  = AN;                   % update AN
 % update residual errors 
   r_c  = r_c  - alpha * Kwc;            Norm = norm(r_c(:))^2;
  for n = 1:nDIM  
   r_drho{n}=r_drho{n}-alpha* Kwdrho{n}; Norm = Norm+norm(r_drho{n}(:))^2;    
  end    
   Error = sqrt(Norm / Norm0); fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
 end % CG iterations
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached:  err/norm = ' num2str(Error)]); end
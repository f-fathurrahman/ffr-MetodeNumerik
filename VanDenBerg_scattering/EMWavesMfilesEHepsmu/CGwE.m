function [w_E] = CGwE(E_inc,FFTG,input)
global nDIM;     CHI_eps = input.CHI_eps;
w_E=cell(1,nDIM);  r_E=cell(1,nDIM);   g_E=cell(1,nDIM);  v_E=cell(1,nDIM); 

itmax  = 1000;   Errcri = 1e-3;  it = 0;     % initialization of iteration
Norm0  = 0;
 for n = 1:nDIM; 
   w_E{n} = zeros(size(E_inc{n})); 
   r_E{n} = CHI_eps .* E_inc{n};         Norm0 = Norm0 + norm(r_E{n}(:))^2;
 end;  
CheckAdjointE(r_E,FFTG,input);          % Check once

Error  = 1;    fprintf('Error =         %g',Error);
while (it < itmax) && ( Error > Errcri) && (Norm0 > eps);
% determine conjugate gradient direction v 
  dummy=cell(1,nDIM);  for n=1:nDIM;  dummy{n}=conj(CHI_eps).*r_E{n}; end;
  KE  = AdjKopE(dummy, FFTG,input); 
  AN  = 0;   
  for n = 1:nDIM;
     g_E{n} =  r_E{n} - KE{n};            %  and window for small CHI_eps 
     g_E{n} = (abs(CHI_eps) >= Errcri).*g_E{n};   AN=AN+norm(g_E{n}(:))^2;
  end
  if it == 0; 
     for n = 1:nDIM;  v_E{n} = g_E{n}; end;
  else
     for n = 1:nDIM;  v_E{n} = g_E{n} + AN/AN_1 * v_E{n};  end;
  end    
% determine step length alpha        
  KE =  KopE(v_E,FFTG,input);
  BN = 0;    
  for n = 1:nDIM;   
      KE{n}= v_E{n} - CHI_eps .* KE{n};       BN = BN + norm(KE{n}(:))^2;  
  end;
  alpha = AN / BN;
% update contrast sources and AN 
  for n = 1:nDIM;  
      w_E{n} = w_E{n} + alpha * v_E{n}; 
  end
  AN_1  = AN;
% update residual errors 
   Norm = 0;
   for n = 1:nDIM;  
     r_E{n} = r_E{n} - alpha * KE{n};     Norm = Norm + norm(r_E{n}(:))^2; 
   end
   Error = sqrt(Norm / Norm0);       fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
 end; % CG iterations
  
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached: err/norm = ' num2str(Error)]); end;
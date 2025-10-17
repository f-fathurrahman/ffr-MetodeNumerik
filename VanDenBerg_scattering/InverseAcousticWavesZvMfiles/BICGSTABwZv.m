function [w_Zv]  = BICGSTABwZv(Zv_inc,input)
global nDIM;    
w_Zv=cell(1,nDIM);   r_Zv=cell(1,nDIM);   v_Zv=cell(1,nDIM); 
Data_Zv{1} = input.CHI_rho .* Zv_inc{1}; 
Data_Zv{2} = input.CHI_rho .* Zv_inc{2};
itmax   = 200;   
it      = 0;                           % initialization of iteration     
for n = 1:nDIM;    w_Zv{n} = zeros(size(Zv_inc{n}));  end                      
for n = 1:nDIM;    r_Zv{n} = Data_Zv{n};              end
Norm_D = norm(r_Zv{1}(:))^2 + norm(r_Zv{2}(:))^2;           
eta_D   = 1 / Norm_D;                    % normalization factor 
Error   = 1;                             % error norm initial error 
fprintf('Error =         %g',Error); 
while (it < itmax) && (Error > input.Errcri)
 % determine gradient directions
   AN = sum(r_Zv{1}(:).*conj(Data_Zv{1}(:))) + ...
        sum(r_Zv{2}(:).*conj(Data_Zv{2}(:)));               
    if it == 0    
       for n = 1:nDIM;  v_Zv{n} = r_Zv{n};  end
    else       
       for n = 1:nDIM;  v_Zv{n} = r_Zv{n} + (AN/AN_1) * v_Zv{n}; end
    end
    AN_1 = AN;
 % determine step length alpha and update residual error r_D
   [Kv_Zv] = KOPwZv(v_Zv,input);              
   for n = 1:nDIM; Kv_Zv{n} = v_Zv{n} - input.CHI_rho .* Kv_Zv{n}; end
         BN  = sum(Kv_Zv{1}(:).*conj(Data_Zv{1}(:))) + ...
               sum(Kv_Zv{2}(:).*conj(Data_Zv{2}(:)));
       alpha = AN / BN; 
   for n = 1:nDIM;  r_Zv{n} = r_Zv{n} - alpha * Kv_Zv{n};  end  
   % + successive overrelaxation (first step of GMRES)
   [Kr_Zv] = KOPwZv(r_Zv,input);  
   for n=1:nDIM; Kr_Zv{n} = r_Zv{n} - input.CHI_rho .* Kr_Zv{n};  end
   beta = ( sum(r_Zv{1}(:).*conj(Kr_Zv{1}(:)))  + ...
             sum(r_Zv{2}(:).*conj(Kr_Zv{2}(:)))  ) ...
            /(norm(Kr_Zv{1}(:))^2 + norm(Kr_Zv{2}(:))^2); 
   % update contrast source w ,v and the residual error r_Zv
   for n=1:nDIM
      w_Zv{n} = w_Zv{n} + alpha * v_Zv{n} + beta * r_Zv{n};  
      v_Zv{n} = (alpha / beta) * (v_Zv{n} - beta * Kv_Zv{n});
      r_Zv{n} = r_Zv{n} - beta * Kr_Zv{n};
   end
   Norm_D = norm(r_Zv{1}(:))^2 + norm(r_Zv{2}(:))^2 ; 
   Error  = sqrt(eta_D * Norm_D); fprintf('\b\b\b\b\b\b\b\b%6f',Error);
   it = it+1; 
end % while
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax  
   disp(['itmax was reached:   err/norm = ' num2str(Error)]); 
end 
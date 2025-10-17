function [dw]  = BICGSTABdw(Pinc,input)
global nDIM;
dw = cell(1,nDIM); dr = cell(1,nDIM); dv = cell(1,nDIM);
Data_dw{1} = input.Rfl{1} .* Pinc{1}; 
Data_dw{2} = input.Rfl{2} .* Pinc{2};
itmax   = 200;   
it      = 0;                             % initialization of iteration
for n = 1:nDIM;  dw{n} = zeros(size(Pinc{n})); end
for n = 1:nDIM;  dr{n} = Data_dw{n};           end  
Norm_D  = norm(dr{1}(:))^2 + norm(dr{2}(:))^2;           
eta_D   = 1 / Norm_D;                    % normalization factor 
Error   = 1;                             % error norm initial error 
fprintf('Error =         %g',Error); 
while (it <= itmax) && (Error > input.Errcri)
 % Determine gradient directions
   AN = sum(dr{1}(:).*conj(Data_dw{1}(:))) + ...
        sum(dr{2}(:).*conj(Data_dw{2}(:)));               
      if it == 0 
         for n = 1:nDIM; dv{n} = dr{n};  end 
      else       
         for n = 1:nDIM; dv{n} = dr{n} +  (AN/AN_1) * dv{n}; end
      end
   AN_1  = AN;
 % determine step length alpha and update residual error r_D
   [Kdv] = KOPdw(dv,input);               
   for n = 1:nDIM; Kdv{n}  = dv{n} - input.Rfl{n} .* Kdv{n};  end
        BN  = sum(Kdv{1}(:).*conj(Data_dw{1}(:))) + ...
              sum(Kdv{2}(:).*conj(Data_dw{2}(:)));
      alpha = AN / BN; 
     for n = 1:nDIM; dr{n}  = dr{n} - alpha * Kdv{n};  end  
    % + successive overrelaxation (first step of GMRES)
    [Kdr] = KOPdw(dr,input); 
    for n = 1:nDIM; Kdr{n} = dr{n} - input.Rfl{n} .* Kdr{n};  end
    beta   = ( sum(dr{1}(:).*conj(Kdr{1}(:)))  + ...
               sum(dr{2}(:).*conj(Kdr{2}(:)))  ) ...
              /(norm(Kdr{1}(:))^2 + norm(Kdr{2}(:))^2); 
    % update contrast source w ,v and the residual error r_D
     for n = 1:nDIM
        dw{n} = dw{n}  + alpha *  dv{n} + beta *  dr{n};   
        dv{n} = (alpha / beta) * (dv{n} - beta * Kdv{n});  
        dr{n}  = dr{n} - beta * Kdr{n};
     end
     Norm_D = norm(dr{1}(:))^2 + norm(dr{2}(:))^2 ; 
     Error  = sqrt(eta_D * Norm_D);
     fprintf('\b\b\b\b\b\b\b\b%6f',Error);       it = it+1; 
end % while
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax  
   disp(['itmax was reached:   err/norm = ' num2str(Error)]); 
end 
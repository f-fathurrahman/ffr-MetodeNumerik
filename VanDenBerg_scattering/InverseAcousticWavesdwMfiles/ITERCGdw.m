function [dw] = ITERCGdw(Pinc,input) 
global nDIM;      Rfl = input.Rfl; 

dw = cell(1,nDIM); dr = cell(1,nDIM); dg = cell(1,nDIM); dv = cell(1,nDIM); 
itmax = 200;  Errcri = input.Errcri;  it = 0;  % initialization iteration                                         
Norm_D = 0;
for n = 1:nDIM 
   dw{n} = zeros(size(Pinc{n}));  
   dr{n} = Rfl{n}.*Pinc{n};            Norm_D = Norm_D + norm(dr{n}(:))^2; 
end   
eta_D = 1 /Norm_D;
Error = 1;  fprintf('Error =         %g',Error);

while (it < itmax) && ( Error > Errcri) 
 % determine conjugate gradient direction v
   dummy=cell(1,nDIM); for n=1:nDIM; dummy{n} = conj(Rfl{n}).*dr{n}; end
   [Kdf] = AdjKOPdw(dummy,input); 
   AN = 0;
   for n = 1:nDIM  
    dg{n} = dr{n} - Kdf{n}; % window for negligible Rfl !!
    dg{n} = (abs(Rfl{n}) >= Errcri) .* dg{n};  
    AN = AN + norm(dg{n}(:))^2;
   end 
   if it == 0
      for n=1:nDIM;  dv{n} = dg{n}; end
   else
      for n=1:nDIM;  dv{n} = dg{n} + AN/AN_1 * dv{n}; end
   end    
% determine step length alpha        
  [Kdv] = KOPdw(dv,input); 
  BN = 0; 
  for n = 1:nDIM   
      Kdv{n} = dv{n} - Rfl{n} .* Kdv{n};     
      BN = BN + norm(Kdv{n}(:))^2;                                   
  end
   alpha = AN / BN;
 % update contrast sources, AN  and update residual errors 
   Norm_D = 0;
   for n = 1:nDIM  
      dw{n} = dw{n} + alpha * dv{n};    
      dr{n} = dr{n} - alpha * Kdv{n};  
      Norm_D = Norm_D + norm(dr{n}(:))^2; 
   end  
   AN_1  = AN;    
   Error = sqrt(eta_D * Norm_D); fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
end % CG iterations
 
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax; disp(['itmax reached:  err/norm = ' num2str(Error)]); end
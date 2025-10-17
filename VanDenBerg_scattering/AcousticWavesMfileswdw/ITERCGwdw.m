function [w,dw] = ITERCGwdw(p_inc,Pinc,input) 
global nDIM;      CHI = input.CHI;   Rfl = input.Rfl; 
dw = cell(1,nDIM); dr = cell(1,nDIM); dg = cell(1,nDIM); dv = cell(1,nDIM); 

itmax = 1000;  Errcri = input.Errcri;  it = 0;  % initialization iteration                                   
    w = zeros(size(p_inc));          
    r = CHI.*p_inc;                      Norm0 = norm(r(:))^2;
for n = 1:nDIM 
   dw{n} = zeros(size(Pinc{n}));  
   dr{n} = Rfl{n}.*Pinc{n};              Norm0 = Norm0 + norm(dr{n}(:))^2; 
end    
CheckAdjointwdw(r,dr,input);     % Check once

Error = 1;   fprintf('Error =         %g',Error);
 while (it < itmax) && ( Error > Errcri) && (Norm0 > eps)
 % determine conjugate gradient direction v
   dummy=cell(1,nDIM); for n = 1:nDIM; dummy{n} = conj(Rfl{n}).*dr{n}; end
   [Kf,Kdf] = AdjKOPwdw(conj(CHI).*r, dummy,input); 
          g =  r - Kf;     % window for negligible CHI !!
          g = (abs(CHI) >= Errcri) .* g;      AN = norm(g(:))^2;
   for n = 1:nDIM  
    dg{n}= dr{n} - Kdf{n}; % window for negligible Rfl !!
    dg{n}= (abs(Rfl{n}) >= Errcri) .* dg{n};  AN = AN + norm(dg{n}(:))^2;
   end 
   if it == 0;      v   =  g;    else  v   =  g    + AN/AN_1 *  v;    end
   for n=1:nDIM   
     if it == 0;  dv{n} = dg{n}; else dv{n}= dg{n} + AN/AN_1 * dv{n}; end
   end    
% determine step length alpha        
  [Kv,Kdv] = KOPwdw(v,dv,input); 
   Kv = v - CHI .* Kv;                        BN = norm(Kv(:))^2; 
   for n = 1:nDIM   
      Kdv{n} = dv{n} - Rfl{n} .* Kdv{n};      BN = BN + norm(Kdv{n}(:))^2;                                   
   end
   alpha = AN / BN;
 % update contrast sources, AN  and update residual errors 
   w  = w + alpha *  v; 
   r  = r - alpha * Kv;                  Norm = norm(r(:))^2;
   for n = 1:nDIM  
      dw{n} = dw{n} + alpha * dv{n};    
      dr{n} = dr{n} - alpha * Kdv{n};    Norm = Norm + norm(dr{n}(:))^2; 
   end  
   AN_1  = AN;    
   Error = sqrt(Norm / Norm0); fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
 end % CG iterations
 
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached:  err/norm = ' num2str(Error)]); end;
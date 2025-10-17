function [w,dw] = ITERCGEwdw(E_inc,dEinc,input)
global nDIM;       CHI = input.CHI_eps;                  Rfl = input.Rfl;
w=cell(1,nDIM);  r=cell(1,nDIM);  dw=cell(1,nDIM);  dr=cell(1,nDIM);  
g=cell(1,nDIM);  v=cell(1,nDIM);  dg=cell(1,nDIM);  dv=cell(1,nDIM);

itmax = 1000; Errcri = input.Errcri; it = 0; Norm0  = 0;  % initialization 
for n = 1:nDIM; 
   w{n} = zeros(size(E_inc{n}));  
  dw{n} = zeros(size(dEinc{n})); 
   r{n} =  CHI   .*  E_inc{n};        Norm0 = Norm0 + norm( r{n}(:))^2; 
  dr{n} = Rfl{n} .* dEinc{n};         Norm0 = Norm0 + norm(dr{n}(:))^2;
end;  
CheckAdjointEwdw(r,dr,input);                       % Check once

Error  = 1;    fprintf('Error =         %g',Error);
while (it < itmax) && ( Error > Errcri) && (Norm0 > eps);
% determine conjugate gradient direction v 
  dummy=cell(1,nDIM); for n=1:nDIM;  dummy{n} = conj(CHI)   .* r{n}; end;
 ddummy=cell(1,nDIM); for n=1:nDIM; ddummy{n} = conj(Rfl{n}).*dr{n}; end;
  [Kf,Kdf]  = AdjKopEwdw(dummy,ddummy,input); AN = 0;
  for n = 1:nDIM;
     g{n} =   r{n} -  Kf{n};    % window for negligible CHI 
     g{n} = (abs(CHI)>=Errcri) .*  g{n};      AN = AN + norm( g{n}(:))^2;
    dg{n} =  dr{n} - Kdf{n};    % window for negligible CHI 
    dg{n} = (abs(CHI)>=Errcri) .* dg{n};      AN = AN + norm(dg{n}(:))^2;
  end;
  for n = 1:nDIM;  
    if it == 0;  v{n} =  g{n}; else  v{n} =  g{n} + AN/AN_1 *  v{n}; end;  
    if it == 0; dv{n} = dg{n}; else dv{n} = dg{n} + AN/AN_1 * dv{n}; end;
  end    
% determine step length alpha        
  [Kv,Kdv] =  KopEwdw(v,dv,input);        BN = 0;    
  for n = 1:nDIM;   
      Kv{n} = v{n} - CHI    .*  Kv{n};    BN = BN + norm( Kv{n}(:))^2;
     Kdv{n} = dv{n} - Rfl{n} .* Kdv{n};   BN = BN + norm(Kdv{n}(:))^2;
  end;
  alpha = AN / BN;                        Norm = 0;
% update contrast sources, AN and resdual errors 
  for n = 1:nDIM;  
      w{n} =  w{n} + alpha *  v{n}; 
     dw{n} = dw{n} + alpha * dv{n}; 
      r{n} =  r{n} - alpha *  Kv{n};      Norm = Norm + norm( r{n}(:))^2;
     dr{n} = dr{n} - alpha * Kdv{n};      Norm = Norm + norm(dr{n}(:))^2; 
  end;
   AN_1 = AN;
  Error = sqrt(Norm / Norm0);       fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
     it = it+1;
 end; % CG iterations
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached: err/norm = ' num2str(Error)]); end;
function [w] = BiCGSTABw(u_inc,input)
% CG_FFT scheme for contrast source integral equation 
CHI = input.CHI;  FFTG = input.FFTG;
 
itmax  = 200;   
Errcri = input.Errcri;  
it     = 0;                              % initialization of iteration
w      = zeros(size(u_inc));             % first guess for contrast source                  
r      = CHI.*u_inc;                     % first residual vector
Norm0  = norm(r(:));                     % normalization factor  
Error  = 1;                              % error norm initial error 
r_tld  = r;

fprintf('Error =         %g',Error);
  
while (it < itmax) && (Error > Errcri) 
    % Determine gradient directions
       AN  = sum(r(:).*conj(r_tld(:)));               
      if it == 0
          v = r; 
      else         
          v = r + (AN/AN_1) * v;
      end
    % determine step length alpha and update residual error r
      Kv    = v - CHI.* Kop(v,FFTG);
      BN    = sum(Kv(:).*conj(r_tld(:)));
      alpha = AN / BN;
      r     = r - alpha * Kv;    
    % + successive overrelaxation (first step of GMRES)
      Kr    = r - CHI.* Kop(r,FFTG);
      beta  = sum(r(:).*conj(Kr(:))) / norm(Kr(:))^2;
    % update contrast source w and AN 
      w     = w + alpha * v + beta * r;
      v     = (alpha/beta) * ( v - beta * Kv);
      AN_1  = AN;  
    % update the residual error r
      r     = r - beta * Kr; 
      Error = norm(r(:)) / Norm0; 
      fprintf('\b\b\b\b\b\b\b\b%6f',Error);       it = it+1;
end % while
  
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax  
   disp(['itmax was reached:   err/norm = ' num2str(Error)]); 
end 

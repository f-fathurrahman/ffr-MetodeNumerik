function [u] = ITERCGu(u_inc,input)
% CGFFT scheme for wave field integral equation
CHI = input.CHI;  FFTG = input.FFTG;

itmax  = 200;   
Errcri = input.Errcri;              
it     = 0;                                % initialization of iteration   
u      = zeros(size(u_inc));               % first guess for field   
r      = u_inc;                            % first residual vector
Norm0  = norm(r(:));                       % normalization factor
Error  = 1;                                % error norm initial error   
 
% check adjoint operator [I-conj(CHI)K*] via inner product ----------------
         dummy =  r - conj(CHI).*AdjKop(r,FFTG);
       Result1 = sum(r(:).*conj(dummy(:)));
         dummy =  r - Kop(CHI.*r,FFTG);
       Result2 = sum(dummy(:).*conj(r(:)));
       fprintf('Check adjoint: %e\n',abs(Result1-Result2));
% -------------------------------------------------------------------------

fprintf('Error =         %g',Error);
 
while (it <itmax) && (Error > Errcri) 
    % determine conjugate gradient direction v
         g  = r - conj(CHI).*AdjKop(r,FFTG);
        AN  = norm(g(:))^2;
      if it == 0
          v = g; 
      else
          v = g + (AN/AN_1) * v;
      end
    % determine step length alpha
      Kv   = v - Kop(CHI.*v,FFTG);
      BN    = norm(Kv(:))^2;
      alpha = AN / BN;
    % update field and AN
      u     = u + alpha * v;
      AN_1  = AN;  
    % update the residual error r
      r     = r - alpha * Kv;
      Error = sqrt(norm(r(:)) / Norm0);
      fprintf('\b\b\b\b\b\b\b\b%6f',Error);              
      it = it+1;
end % while
  
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax 
   disp(['itmax was reached:   Error = ' num2str(Error)]); 
end

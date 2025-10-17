function [w] = ITERCGw(u_inc,input)
% CG_FFT scheme for contrast source integral equation 
CHI = input.CHI;  FFTG = input.FFTG;
 
itmax  = 200;   
Errcri = input.Errcri;  
it     = 0;                                % initialization of iteration
w      = zeros(size(u_inc));               % first guess for contrast source                  
r      = CHI.*u_inc;                       % first residual vector
Norm0  = norm(r(:));                       % normalization factor  
Error  = 1;                                % error norm initial error 

% check adjoint operator [I-K*conj(CHI)] via inner product ----------------
         dummy = r - AdjKop(conj(CHI).*r,FFTG);
       Result1 = sum(r(:).*conj(dummy(:)));
         dummy = r - CHI.* Kop(r,FFTG);
       Result2 = sum(dummy(:).*conj(r(:)));
       fprintf('Check adjoint: %e\n',abs(Result1-Result2));
% -------------------------------------------------------------------------

fprintf('Error =         %g',Error);
  
while (it < itmax) && (Error > Errcri) 
    % determine conjugate gradient direction v
         g = r - AdjKop(conj(CHI).*r,FFTG);
         g = (abs(CHI) >= Errcri) .* g; % window for negligible contrast!!!
       AN  = norm(g(:))^2;
      if it == 0
          v = g; 
      else
          v = g + (AN/AN_1) * v;
      end
    % determine step length alpha
      Kv    = v - CHI.* Kop(v,FFTG);
      BN    = norm(Kv(:))^2;
      alpha = AN / BN;
    % update contrast source w and AN 
      w     = w + alpha * v; 
      AN_1  = AN;  
    % update the residual error r
      r     = r - alpha * Kv; 
      Error =sqrt(norm(r(:)) / Norm0);
      fprintf('\b\b\b\b\b\b\b\b%6f',Error);       it = it+1;
end % while
  
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax
   disp(['itmax was reached:   err/norm = ' num2str(Error)]); 
end 

function [w] = ITERCGw(u_inc,input)
% CG_FFT scheme for contrast source integral equation 
CHI = input.CHI;  FFTG = input.FFTG;
 
itmax   = 200;   
Errcri  = input.Errcri;  
it      = 0;                              % initialization of iteration
w       = zeros(size(u_inc));             % first guess for contrast source                  
r_D     = CHI.*u_inc;                     % first residual vector
Norm_D  = norm(r_D(:))^2;                    
eta_D   = 1 / Norm_D;                     % normalization factor
Error   = 1;                              % error norm initial error 
fprintf('Error =         %g',Error);

% check adjoint operator [I-K*conj(CHI)] via inner product ----------------
%          dummy = r_D - AdjKop(conj(CHI).*r_D,FFTG);
%        Result1 = sum(r_D(:).*conj(dummy(:)));
%          dummy = r_D - CHI.* Kop(r_D,FFTG);
%        Result2 = sum(dummy(:).*conj(r_D(:)));
%        fprintf('Check adjoint: %e\n',abs(Result1-Result2));
% -------------------------------------------------------------------------

while (it < itmax) && (Error > Errcri) 
    % determine conjugate gradient direction v
         g = r_D - AdjKop(conj(CHI).*r_D,FFTG);
         g = (abs(CHI) >= Errcri) .* g; % window for negligible contrast!!
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
      r_D    = r_D - alpha * Kv; 
      Norm_D = norm(g(:))^2;
      Error  = sqrt(eta_D * Norm_D);
      fprintf('\b\b\b\b\b\b\b\b%6f',Error);                it = it+1;
end % while
  
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of iterations is ' num2str(it)]);
if it == itmax  
   disp(['itmax was reached:   err/norm = ' num2str(Error)]); 
end 
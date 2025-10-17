function [w] = BICGSTABw(u_inc,input)
CHI = input.CHI;   Data_D = CHI.*u_inc;    FFTG = input.FFTG;

itmax  = 200;     
Errcri = input.Errcri;  
it     = 0;                              % initialization of iteration
w      = zeros(size(u_inc));             % first guess for contrast source                  
r_D    = Data_D;                         % first residual vector
Norm_D = norm(r_D(:));                   % normalization factor  
eta_D  = 1 / Norm_D;
Error  = 1;                              % error norm initial error 
  
while (it < itmax) && (Error > Errcri)
    
    % Determine gradient directions
       AN  = sum(r_D(:).*conj(Data_D(:)));               
      if it == 0
          v = r_D; 
      else       
          v = r_D +  (AN/AN_1) * v;
      end
    % determine step length alpha and update residual error r_D
      Kv    = v - CHI.* Kop(v,FFTG);
      BN    = sum(Kv(:).*conj(Data_D(:)));
      alpha = AN / BN;
      r_D   = r_D - alpha * Kv;    
    % + successive overrelaxation (first step of GMRES)
      Kr    = r_D - CHI.* Kop(r_D,FFTG);
      beta  = sum(r_D(:).*conj(Kr(:))) / norm(Kr(:))^2; 
    % update contrast source w ,v and AN 
      w     = w + alpha * v   + beta * r_D;  
      v     = (alpha / beta) * (v - beta * Kv);
      AN_1  = AN;  
    % update the residual error r_D
      r_D   = r_D - beta * Kr;
      Norm_D = norm(r_D(:))^2;
      Error = sqrt(eta_D * Norm_D);     
      it = it+1;
end % while
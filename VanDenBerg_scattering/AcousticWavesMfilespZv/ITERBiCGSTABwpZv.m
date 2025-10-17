function [w_p,w_Zv] = ITERBiCGSTABwpZv(p_inc,Zv_inc,input)
% GMRES_FFT scheme for contrast source integral equation Aw = b
  global nDIM;     [N,~] = size(input.CHI_kap(:)); 
  itmax  = 1000;  Errcri = input.Errcri;  
                                                         % known vector
                 b(1:N,1)       = input.CHI_kap(:) .* p_inc(:); 
                 b(N+1:2*N,1)   = input.CHI_rho(:) .* Zv_inc{1}(:); 
  if nDIM >= 2;  b(2*N+1:3*N,1) = input.CHI_rho(:) .* Zv_inc{2}(:); end
  if nDIM == 3;  b(3*N+1:4*N,1) = input.CHI_rho(:) .* Zv_inc{3}(:); end
  
  w = bicgstab(@(w) Aw(w,input), b, Errcri, itmax);      % call BICGSTAB
 
  [w_p,w_Zv] = vector2matrix(w,input);                   % output matrices
end %----------------------------------------------------------------------


function y = Aw(w,input) 
  global nDIM;      [N,~] = size(input.CHI_kap(:)); 

  [w_p,w_Zv] = vector2matrix(w,input); 
  [Kp,KZv]   = KOPpZv(w_p,w_Zv,input); 
          Kp = w_p - input.CHI_kap .* Kp;
  for n = 1:nDIM   
      KZv{n} = w_Zv{n} - input.CHI_rho .* KZv{n};  
  end
                 y(1:N,1)       = Kp(:);
                 y(N+1:2*N,1)   = KZv{1}(:);
  if nDIM >= 2;  y(2*N+1:3*N,1) = KZv{2}(:);  end
  if nDIM == 3;  y(3*N+1:4*N,1) = KZv{3}(:);  end
end %----------------------------------------------------------------------


function [w_p,w_Zv] = vector2matrix(w,input)
% Modify vector output from 'bicgstab' to matrices for further computation
  global nDIM;      [N,~] = size(input.CHI_kap(:));   w_Zv = cell(1,nDIM); 

  if nDIM == 2;  DIM = [input.N1,input.N2];           end
  if nDIM == 3;  DIM = [input.N1,input.N2,input.N3];  end

                     w_p = w(1:N,1); 
                 w_Zv{1} = w(N+1:2*N,1);
  if nDIM >= 2;  w_Zv{2} = w(2*N+1:3*N,1);  end
  if nDIM == 3;  w_Zv{3} = w(3*N+1:4*N,1);  end
 
  if nDIM >= 2;     w_p = reshape(w_p,DIM); end
  for n = 1:nDIM
     if nDIM >= 2;  w_Zv{n} = reshape(w_Zv{n},DIM);  end          
  end
  
end
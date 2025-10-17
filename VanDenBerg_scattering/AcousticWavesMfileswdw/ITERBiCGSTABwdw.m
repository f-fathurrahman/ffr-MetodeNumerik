function [wc,dw] = ITERBiCGSTABwdw(p_inc,Pinc,input)
  % GMRES_FFT scheme for contrast source integral equation Aw = b
  global nDIM;       [N,~] = size(input.CHI(:)); 
  itmax  = 1000;    Errcri = input.Errcri; 
                                                      % known vector
                 b(1:N,1)       = input.CHI(:)    .* p_inc(:); 
                 b(N+1:2*N,1)   = input.Rfl{1}(:) .* Pinc{1}(:); 
  if nDIM >= 2;  b(2*N+1:3*N,1) = input.Rfl{2}(:) .* Pinc{2}(:); end
  if nDIM == 3;  b(3*N+1:4*N,1) = input.Rfl{3}(:) .* Pinc{3}(:); end
 
w = bicgstab(@(w) Aw(w,input), b, Errcri, itmax);     % call BICGSTAB

[wc,dw] = vector2matrix(w,input);                     % output matrices    
end  %---------------------------------------------------------------------


function y = Aw(w,input) 
  global nDIM;      [N,~] = size(input.CHI(:)); 

  [wc,dw]  = vector2matrix(w,input);    
  [Kw,Kdw] = KOPwdw(wc,dw,input); 
        Kw = wc - input.CHI .* Kw;
  for n = 1:nDIM   
    Kdw{n} = dw{n} - input.Rfl{n} .* Kdw{n};  
  end
                 y(1:N,1)       = Kw(:);
                 y(N+1:2*N,1)   = Kdw{1}(:);
  if nDIM >= 2;  y(2*N+1:3*N,1) = Kdw{2}(:);  end
  if nDIM == 3;  y(3*N+1:4*N,1) = Kdw{3}(:);  end
end %----------------------------------------------------------------------


function [wc,dw] = vector2matrix(w,input)
% Modify vector output from 'bicgstab' to matrices for further computation
  global nDIM;      [N,~] = size(input.CHI(:));   dw = cell(1,nDIM); 

  if nDIM == 2;    DIM = [input.N1,input.N2];           end
  if nDIM == 3;    DIM = [input.N1,input.N2,input.N3];  end

                    wc = w(1:N,1); 
                 dw{1} = w(N+1:2*N,1);
  if nDIM >= 2;  dw{2} = w(2*N+1:3*N,1);  end
  if nDIM == 3;  dw{3} = w(3*N+1:4*N,1);  end
 
  if nDIM >= 2;     wc = reshape(wc,DIM); end
  for n = 1:nDIM
     if nDIM >= 2; dw{n} = reshape(dw{n},DIM);   end 
  end  
  
end 
function [w_E] = ITERBiCGSTABwE(E_inc,input)
% BiCGSTAB scheme for contrast source integral equation Aw = b
  global nDIM;     
  itmax  = 1000;  Errcri = input.Errcri;  [N,~] = size(input.CHI_eps(:)); 

          b(    1:  N,1) = input.CHI_eps(:) .* E_inc{1}(:);        
          b(  N+1:2*N,1) = input.CHI_eps(:) .* E_inc{2}(:);
          
  if nDIM == 3  
          b(2*N+1:3*N,1) = input.CHI_eps(:) .* E_inc{3}(:); 
  end % if                                   
      
   w = bicgstab(@(w) Aw(w,input), b, Errcri, itmax);     % call BICGSTAB

  [w_E] = vector2matrix(w,input);                        % output matrices
 
end %---------------------------------------------------------------------- 


function y = Aw(w,input) 
  global nDIM;      [N,~] = size(input.CHI_eps(:));  
  
  [w_E]  = vector2matrix(w,input);     
  [Kw_E] = KopE(w_E,input); 

          y(    1:  N,1) = w_E{1}(:) - input.CHI_eps(:) .* Kw_E{1}(:);  
          y(  N+1:2*N,1) = w_E{2}(:) - input.CHI_eps(:) .* Kw_E{2}(:);
          
  if nDIM == 3  
          y(2*N+1:3*N,1) = w_E{3}(:) - input.CHI_eps(:) .* Kw_E{3}(:);    
  end % if  
  
end %----------------------------------------------------------------------


function [w_E] = vector2matrix(w,input)
% Modify vector output from 'bicgstab' to matrices for further computation
  global nDIM;   [N,~] = size(input.CHI_eps(:));      w_E = cell(1,nDIM); 
 
  if nDIM == 2;   DIM = [input.N1,input.N2];           end % if
  if nDIM == 3;   DIM = [input.N1,input.N2,input.N3];  end % if

             w_E{1} = reshape(w(    1:  N,1),DIM); 
             w_E{2} = reshape(w(  N+1:2*N,1),DIM);  
             
  if nDIM == 3  
             w_E{3} = reshape(w(2*N+1:3*N,1),DIM);   
  end %if
 
end
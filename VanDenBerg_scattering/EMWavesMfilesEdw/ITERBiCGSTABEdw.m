function [dw] = ITERBiCGSTABEdw(E_inc,input)
% BiCGSTAB scheme for contrast source integral equation Aw = b
  global nDIM; 
  itmax = 1000; Errcri = input.Errcri; [N,~] = size(input.Rfl{1}(:)); 

  if nDIM == 2;       
                      b(    1:  N,1) = input.Rfl{1}(:) .* E_inc{2,1}(:); 
                      b(  N+1:2*N,1) = input.Rfl{2}(:) .* E_inc{1,2}(:);  
  
  elseif nDIM == 3; 
                      b(    1:  N,1) = input.Rfl{1}(:) .* E_inc{2,1}(:); 
                      b(  N+1:2*N,1) = input.Rfl{1}(:) .* E_inc{3,1}(:); 
                      b(2*N+1:3*N,1) = input.Rfl{2}(:) .* E_inc{1,2}(:);  
                      b(3*N+1:4*N,1) = input.Rfl{2}(:) .* E_inc{3,2}(:);
                      b(4*N+1:5*N,1) = input.Rfl{3}(:) .* E_inc{1,3}(:);  
                      b(5*N+1:6*N,1) = input.Rfl{3}(:) .* E_inc{2,3}(:);                   
  end; % if 
 
    w  = bicgstab(@(w) Aw(w,input), b, Errcri, itmax);     % call BICGSTAB
  
  [dw] = Edwvector2matrix(w,input); 
  
end %---------------------------------------------------------------------- 

function y = Aw(w,input) 
  global nDIM;      [N,~] = size(input.CHI_eps(:));  

   [dw] = Edwvector2matrix(w,input);    
  
  [Kdw] = KopEdw(dw,input); 
  
  if nDIM == 2;          
            y(    1:  N,1) = dw{2,1}(:) - input.Rfl{1}(:) .* Kdw{2,1}(:);
            y(  N+1:2*N,1) = dw{1,2}(:) - input.Rfl{2}(:) .* Kdw{1,2}(:); 
          
  elseif nDIM == 3;               
            y(    1:  N,1) = dw{2,1}(:) - input.Rfl{1}(:) .* Kdw{2,1}(:);
            y(  N+1:2*N,1) = dw{3,1}(:) - input.Rfl{1}(:) .* Kdw{3,1}(:); 
            y(2*N+1:3*N,1) = dw{1,2}(:) - input.Rfl{2}(:) .* Kdw{1,2}(:);  
            y(3*N+1:4*N,1) = dw{3,2}(:) - input.Rfl{2}(:) .* Kdw{3,2}(:);
            y(4*N+1:5*N,1) = dw{1,3}(:) - input.Rfl{3}(:) .* Kdw{1,3}(:);  
            y(5*N+1:6*N,1) = dw{2,3}(:) - input.Rfl{3}(:) .* Kdw{2,3}(:);            
  end; %if
end
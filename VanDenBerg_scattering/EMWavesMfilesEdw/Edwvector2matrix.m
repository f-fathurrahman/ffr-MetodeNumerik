function [dw] = Edwvector2matrix(w,input)
% Modify vector output from 'bicgstab' to matrices for further computation
  global nDIM;  
  
  [N,~] = size(input.Rfl{1}(:));    
     dw = cell(nDIM,nDIM); 
  
  if nDIM == 2;                     DIM = [input.N1,input.N2];  

   dw{2,1} = reshape(w(    1:  N,1),DIM);  
   dw{1,2} = reshape(w(  N+1:2*N,1),DIM);

  elseif  nDIM == 3;                DIM = [input.N1,input.N2,input.N3]; 
    
   dw{2,1} = reshape(w(    1:  N,1),DIM);  
   dw{3,1} = reshape(w(  N+1:2*N,1),DIM); 
   dw{1,2} = reshape(w(2*N+1:3*N,1),DIM);
   dw{3,2} = reshape(w(3*N+1:4*N,1),DIM);   
   dw{1,3} = reshape(w(4*N+1:5*N,1),DIM);  
   dw{2,3} = reshape(w(5*N+1:6*N,1),DIM); 
 
  end;  %if  
end
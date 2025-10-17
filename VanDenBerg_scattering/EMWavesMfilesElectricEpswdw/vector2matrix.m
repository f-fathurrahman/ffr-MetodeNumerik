function [w_E,dwE] = vector2matrix(w,input)
% Modify vector output from 'bicgstab' to matrices for further computation
  global nDIM;  
  
  [N,~] = size(input.CHI_eps(:));      
    w_E = cell(1,nDIM);  
    dwE = cell(1,nDIM); 
  
  if nDIM == 2;                         DIM = [input.N1,input.N2];    
     
        w_E{1} = reshape(w(    1:  N,1),DIM); 
        w_E{2} = reshape(w(  N+1:2*N,1),DIM); 
        dwE{1} = reshape(w(2*N+1:3*N,1),DIM);
        dwE{2} = reshape(w(3*N+1:4*N,1),DIM); 
        
  end; %if
  
  if nDIM == 3;                       DIM = [input.N1,input.N2,input.N3];
      
      w_E{1} = reshape(w(    1:  N,1),DIM); 
      w_E{2} = reshape(w(  N+1:2*N,1),DIM); 
      w_E{3} = reshape(w(2*N+1:3*N,1),DIM);
      dwE{1} = reshape(w(3*N+1:4*N,1),DIM);
      dwE{2} = reshape(w(4*N+1:5*N,1),DIM);
      dwE{3} = reshape(w(5*N+1:6*N,1),DIM);  
      
  end; %if
 
end
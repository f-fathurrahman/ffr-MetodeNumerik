function [w_E,w_H] = EHvector2matrix(w,input)
% Modify vector output from 'bicgstab' to matrices for further computation
global nDIM;
[N,~] = size(input.CHI_eps(:));    w_E = cell(1,nDIM); w_H = cell(1,nDIM); 

if nDIM == 2;                      DIM = [input.N1,input.N2];  

   w_E{1} = reshape(w(    1:  N,1),DIM);  
   w_E{2} = reshape(w(  N+1:2*N,1),DIM); 
   w_H{3} = reshape(w(2*N+1:3*N,1),DIM);
      
elseif  nDIM == 3;                 DIM = [input.N1,input.N2,input.N3]; 
    
   w_E{1} = reshape(w(    1:  N,1),DIM);
   w_E{2} = reshape(w(  N+1:2*N,1),DIM);
   w_E{3} = reshape(w(2*N+1:3*N,1),DIM);
   w_H{1} = reshape(w(3*N+1:4*N,1),DIM);  
   w_H{2} = reshape(w(4*N+1:5*N,1),DIM);  
   w_H{3} = reshape(w(5*N+1:6*N,1),DIM);     
end      
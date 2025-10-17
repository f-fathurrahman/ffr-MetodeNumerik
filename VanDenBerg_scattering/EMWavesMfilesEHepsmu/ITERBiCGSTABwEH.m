function [w_E,w_H] = ITERBiCGSTABwEH(E_inc,H_inc,input)
% BiCGSTAB scheme for contrast source integral equation Aw = b
 global nDIM;  
  
 itmax  = 1000;  Errcri = input.Errcri;   [N,~] = size(input.CHI_eps(:)); 

 if nDIM == 2;      b(    1:  N,1) = input.CHI_eps(:) .* E_inc{1}(:);        
                    b(  N+1:2*N,1) = input.CHI_eps(:) .* E_inc{2}(:); 
                    b(2*N+1:3*N,1) = input.CHI_mu(:)  .* H_inc{3}(:);
                  
 elseif nDIM == 3;  b(    1:  N,1) = input.CHI_eps(:) .* E_inc{1}(:);        
                    b(  N+1:2*N,1) = input.CHI_eps(:) .* E_inc{2}(:);                                        
                    b(2*N+1:3*N,1) = input.CHI_eps(:) .* E_inc{3}(:);
                    b(3*N+1:4*N,1) = input.CHI_mu(:)  .* H_inc{1}(:); 
                    b(4*N+1:5*N,1) = input.CHI_mu(:)  .* H_inc{2}(:);
                    b(5*N+1:6*N,1) = input.CHI_mu(:)  .* H_inc{3}(:);
 end;  
 
 w = bicgstab(@(w) Aw(w,input), b, Errcri, itmax);       % call BICGSTAB
  
 [w_E,w_H] = EHvector2matrix(w,input);  
  
end %---------------------------------------------------------------------- 

function y = Aw(w,input) 
 global nDIM;      [N,~] = size(input.CHI_eps(:));  

 [w_E,w_H] = EHvector2matrix(w,input);    
  
 [Kw_E,Kw_H] = KopEH(w_E,w_H,input); 
 
 if nDIM == 2;     
     y(    1:  N,1) = w_E{1}(:) - input.CHI_eps(:) .* Kw_E{1}(:);
     y(  N+1:2*N,1) = w_E{2}(:) - input.CHI_eps(:) .* Kw_E{2}(:);
     y(2*N+1:3*N,1) = w_H{3}(:) - input.CHI_mu(:)  .* Kw_H{3}(:); 
               
 elseif nDIM ==3;  
     y(    1:  N,1) = w_E{1}(:) - input.CHI_eps(:) .* Kw_E{1}(:);
     y(  N+1:2*N,1) = w_E{2}(:) - input.CHI_eps(:) .* Kw_E{2}(:);
     y(2*N+1:3*N,1) = w_E{3}(:) - input.CHI_eps(:) .* Kw_E{3}(:); 
     y(3*N+1:4*N,1) = w_H{1}(:) - input.CHI_mu(:)  .* Kw_H{1}(:); 
     y(4*N+1:5*N,1) = w_H{2}(:) - input.CHI_mu(:)  .* Kw_H{2}(:); 
     y(5*N+1:6*N,1) = w_H{3}(:) - input.CHI_mu(:)  .* Kw_H{3}(:); 
 end;
end 
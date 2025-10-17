function [KE] = KopEdw(dw,input)   
global nDIM;

if nDIM == 2;   
factor = 2 / (input.gamma_0^2 * input.dx);  % for num dif Green function 

Kdw{2,1}  = factor * Kop(dw{2,1},input.FFTG);      % input dw_2 on grid(1)            
Kdw{1,2}  = factor * Kop(dw{1,2},input.FFTG);      % input dw_1 on grid(2)              
 
% Compute tangential component (1) on grid(2)              
  KE{1,2} = gradj( Kdw{1,2}-interpolate(Kdw{2,1}, 2,1,input), 2,input); 
        
% Compute tangential component (2) on grid(1)  
  KE{2,1} = gradj( Kdw{2,1}-interpolate(Kdw{1,2}, 1,2,input), 1,input); 
  
elseif nDIM == 3;  
    
  [KE] = KopEdw3D(dw,input);  
  
end % if        
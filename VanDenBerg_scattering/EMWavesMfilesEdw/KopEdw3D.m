function [KE] = KopEdw3D(dw,input)

factor = 2 / (input.gamma_0^2 * input.dx);  % for num dif Green function

Kdw{2,1}  = factor * Kop(dw{2,1},input.FFTG);      % input dw_2 on grid(1)
Kdw{3,1}  = factor * Kop(dw{3,1},input.FFTG);      % input dw_3 on grid(1)

Kdw{1,2}  = factor * Kop(dw{1,2},input.FFTG);      % input dw_1 on grid(2)              
Kdw{3,2}  = factor * Kop(dw{3,2},input.FFTG);      % input dw_3 on grid(2)

Kdw{1,3}  = factor * Kop(dw{1,3},input.FFTG);      % input dw_1 on grid(3)               
Kdw{2,3}  = factor * Kop(dw{2,3},input.FFTG);      % input dw_2 on grid(3)

% Compute tangential field component 1 on grid(2) ------------------------  
  temp = Kdw{1,2} - interpolate(Kdw{2,1},2,1,input);
            KE{1,2} = gradj(temp,2,input);
  temp = interpolate(Kdw{1,3},2,3,input) - interpolate(Kdw{3,1},2,1,input);                   
  KE{1,2} = KE{1,2} + gradj(temp,3,input);
  
% Compute tangential field component 1 on grid(3) ------------------------  
  temp = Kdw{1,3} - interpolate(Kdw{3,1},3,1,input);
            KE{1,3} = gradj(temp,3,input);
  temp = interpolate(Kdw{1,2},3,2,input) - interpolate(Kdw{2,1},3,1,input);                   
  KE{1,3} = KE{1,3} + gradj(temp,2,input); 
  
  
% Compute tangential field component 2 on grid(3) ------------------------  
  temp = Kdw{2,3} - interpolate(Kdw{3,2},3,2,input);
            KE{2,3} = gradj(temp,3,input);
  temp = interpolate(Kdw{2,1},3,1,input) - interpolate(Kdw{1,2},3,2,input);                   
  KE{2,3} = KE{2,3} + gradj(temp,1,input); 
  
% Compute tangential field component 2 on grid(1) ------------------------  
  temp = Kdw{2,1} - interpolate(Kdw{1,2},1,2,input);
            KE{2,1} = gradj(temp,1,input);
  temp = interpolate(Kdw{2,3},1,3,input) - interpolate(Kdw{3,2},1,2,input);                   
  KE{2,1} = KE{2,1} + gradj(temp,3,input);
  
  
% Compute tangential field component 3 on grid(1) ------------------------  
  temp = Kdw{3,1} - interpolate(Kdw{1,3},1,3,input);
            KE{3,1} = gradj(temp,1,input);
  temp = interpolate(Kdw{3,2},1,2,input) - interpolate(Kdw{2,3},1,3,input);                   
  KE{3,1} = KE{3,1} + gradj(temp,2,input); 
  
% Compute tangential field component 3 on grid(2) ------------------------  
  temp = Kdw{3,2} - interpolate(Kdw{2,3},2,3,input);
            KE{3,2} = gradj(temp,2,input);
  temp = interpolate(Kdw{3,1},2,1,input) - interpolate(Kdw{1,3},2,3,input);                   
  KE{3,2} = KE{3,2} + gradj(temp,1,input);
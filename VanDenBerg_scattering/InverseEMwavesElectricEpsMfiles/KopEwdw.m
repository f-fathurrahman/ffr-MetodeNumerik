function [Kw,Kdw] = KopEwdw(w,dw,input)                 

% Correction factor for interface integrals
  factor = 2 / (input.gamma_0^2 * input.dx); 
 
  Kw1  = Kop(w{1},input.FFTG);   
  Kw2  = Kop(w{2},input.FFTG);    
  Kdw1 = factor * Kop(dw{1},input.FFTG);
  Kdw2 = factor * Kop(dw{2},input.FFTG); 

% Interpolate the tangential components of Kdw to grid 0  
  Kw1 = Kw1 + interpolate(gradj(Kdw1,1,input),0,1,input);
  Kw2 = Kw2 + interpolate(gradj(Kdw2,2,input),0,2,input);  
  
% Extrapolate the normal components of Kdw to grid 0    
  Kw{1} = Kw1 + extrapolate(gradj(Kdw2,1,input),0,2,input);
  Kw{2} = Kw2 + extrapolate(gradj(Kdw1,2,input),0,1,input);
                
% Interpolate the tangential components of Kw to the grid 0
  Kdw{1} = interpolate(Kw{1},1,0,input);
  Kdw{2} = interpolate(Kw{2},2,0,input);                
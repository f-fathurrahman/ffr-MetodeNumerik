function [KZv] = AdjKOPwZv(fZv,input)
global nDIM;  
KZv = cell(1,nDIM);   
gam0 = input.gamma_0; 
FFTG = input.FFTG;  

% Acoustic operator with grad div  differential operator
  KZv{1} = AdjKop(fZv{1},FFTG)/conj(gam0^2); 
  KZv{2} = AdjKop(fZv{2},FFTG)/conj(gam0^2);   
   
  KZv = graddiv(KZv,input);  
       
  KZv{1} = KZv{1} + fZv{1}; 
  KZv{2} = KZv{2} + fZv{2}; 
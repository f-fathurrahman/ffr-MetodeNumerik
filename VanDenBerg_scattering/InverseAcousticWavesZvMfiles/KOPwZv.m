function [KZv] = KOPwZv(wZv,input)
global nDIM;
KZv  = cell(1,nDIM);
gam0 = input.gamma_0; 
FFTG = input.FFTG; 

% Acoustic operator with grad div  differential operator
  KZv{1} = Kop(wZv{1},FFTG)/gam0^2; 
  KZv{2} = Kop(wZv{2},FFTG)/gam0^2; 
  
  KZv = graddiv(KZv,input); 
         
  KZv{1} =  KZv{1} + wZv{1}; 
  KZv{2} =  KZv{2} + wZv{2}; 
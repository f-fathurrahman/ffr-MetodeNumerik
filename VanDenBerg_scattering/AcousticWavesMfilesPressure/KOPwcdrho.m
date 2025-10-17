function [Kwdw] = KOPwcdrho(w,dw,input)
global nDIM;   

switch lower(input.method)
    
case 'analytical'         % Acoustic operator using analytical derivatives 
   Kwdw = Kop(w,input.FFTG);
   for n = 1:nDIM  
       Kwdw = Kwdw - Kop(dw{n},input.FFTdG{n}); 
   end

case 'numerical'          % Acoustic operator using numerical derivatives 
   gam0  = input.gamma_0;   Kwdw = cell(1,nDIM);
   Kw = Kop(w,input.FFTG);  
   for n = 1:nDIM    
       Kwdw{n} = Kop(dw{n},input.FFTG); 
   end
   Kwdw = Kw  - DIV(Kwdw,input) / gam0^2;  
   
end % switch 
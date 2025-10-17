function [Kp,KZv] = KOPpZv(wp,wZv,input)
global nDIM;
gam0=input.gamma_0; FFTG=input.FFTG; FFTdG=input.FFTdG;  KZv=cell(1,nDIM);

switch lower(input.method)
case 'analytical'
   % Acoustic operator using analytical derivatives  
   Kp = Kop(wp,FFTG); 
   for n = 1:nDIM
       Kp = Kp - gam0 *Kop(wZv{n},FFTdG{n});
       KZv{n} = Kop(wZv{n},FFTG)/gam0^2;    
   end  
   KZv = graddiv(KZv,input); 
   for n = 1:nDIM        
       KZv{n} =  - gam0 * Kop(wp,FFTdG{n}) + KZv{n} + wZv{n}; 
   end

case 'numerical'
   % Acoustic operator using numerical derivatives 
   Kp = Kop(wp,FFTG); 
   for n=1:nDIM
      KZv{n} = Kop(wZv{n},FFTG);          % use KZv als temporary storage
   end
   Kp = Kp - DIV(KZv,input) / gam0;
   KZv = GRAD(-Kp/gam0,input);
   for n=1:nDIM
       KZv{n} = KZv{n} + wZv{n};
   end

end % switch
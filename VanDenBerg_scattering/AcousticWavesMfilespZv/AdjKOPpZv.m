function [Kp,KZv] = AdjKOPpZv(fp,fZv,input)
global nDIM;     
gam0=input.gamma_0; FFTG=input.FFTG; FFTdG=input.FFTdG;  KZv=cell(1,nDIM);

switch lower(input.method)
case 'analytical'
   % Acoustic operator using analytical derivatives 
   Kp = AdjKop(fp,FFTG); 
   for n = 1:nDIM
       Kp = Kp - conj(gam0) * AdjKop(fZv{n},FFTdG{n});
       KZv{n} = AdjKop(fZv{n},FFTG)/conj(gam0^2);     
   end  
   KZv = graddiv(KZv,input);  
   for n = 1:nDIM        
       KZv{n} = - conj(gam0) * AdjKop(fp,FFTdG{n}) + KZv{n} + fZv{n}; 
   end

case 'numerical'
   % Acoustic operator using numerical derivatives 
   Kp = AdjKop(fp,FFTG);
   for n=1:nDIM
       KZv{n} = AdjKop(fZv{n},FFTG);        % use KZv as temporary storage
   end
   Kp = Kp + DIV(KZv,input) / conj(gam0); 
   KZv = GRAD(Kp/conj(gam0),input); 
   for n=1:nDIM
       KZv{n} = KZv{n} + fZv{n};    
   end

end % switch
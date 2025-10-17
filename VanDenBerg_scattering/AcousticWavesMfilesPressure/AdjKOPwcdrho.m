function [Kf,dKf] = AdjKOPwcdrho(f,df,input)
global nDIM;        dKf = cell(1,nDIM); 

switch lower(input.method)
    
case 'analytical'      % Acoustic operator using analytical derivatives   
   for n = 1:nDIM
       f = f + df{n};
   end
   Kf = AdjKop(f,input.FFTG);
 
for n = 1:nDIM
  dKf{n} = -AdjKop(f,input.FFTdG{n}); 
end

case 'numerical'       % Acoustic operator using numerical derivatives
   gam0 = input.gamma_0; 
   for n = 1:nDIM
       f = f + df{n};
   end
   Kf = AdjKop(f,input.FFTG);    
  dKf = GRAD(Kf/conj(gam0^2),input); 
  
end % switch
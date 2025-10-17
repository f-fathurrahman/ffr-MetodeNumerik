function [Kdw] = KOPdw(dw,input)   
   
FFTdG{1} = input.FFTdG{1} * 2 / input.dx;  
FFTdG{2} = input.FFTdG{2} * 2 / input.dx; 

Kdiag  = Kop(dw{1},FFTdG{1}); % ----------------------------------------                                      
Kdw{1} = Kdiag;                          
Kdw{2} = interpolate(Kdiag,2,1,input);  

Kdiag  = Kop(dw{2},FFTdG{2}); % ----------------------------------------                                     
Kdw{1} = Kdw{1} + interpolate(Kdiag,1,2,input);  
Kdw{2} = Kdiag  + Kdw{2};            
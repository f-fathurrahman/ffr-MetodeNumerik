function [KE,KH] = KopEH(wE,wH,input) 
global nDIM;
gam0  = input.gamma_0; 
FFTG  = input.FFTG;
 
if nDIM == 2                             
     KE{1} = Kop(wE{1},FFTG);               KH{1} = zeros(size(wH{3}));
     KE{2} = Kop(wE{2},FFTG);               KH{2} = zeros(size(wH{3}));
     KE{3} = zeros(size(wH{3}));            KH{3} = Kop(wH{3},FFTG);  
 graddiv_E = graddiv(KE,input);            curl_E = CURL(KE,input); 
                                           curl_H = CURL(KH,input); 
     KE{1} = KE{1} - graddiv_E{1}/gam0^2 + curl_H{1}/gam0;   
     KE{2} = KE{2} - graddiv_E{2}/gam0^2 + curl_H{2}/gam0;
     KH{3} = KH{3} - curl_E{3}/gam0; 
    
elseif nDIM == 3    
     KE{1} = Kop(wE{1},FFTG);               KH{1} = Kop(wH{1},FFTG);
     KE{2} = Kop(wE{2},FFTG);               KH{2} = Kop(wH{2},FFTG);
     KE{3} = Kop(wE{3},FFTG);               KH{3} = Kop(wH{3},FFTG);
 graddiv_E = graddiv(KE,input);            curl_E = CURL(KE,input);
 graddiv_H = graddiv(KH,input);            curl_H = CURL(KH,input);
     KE{1} = KE{1} - graddiv_E{1}/gam0^2 + curl_H{1}/gam0;   
     KE{2} = KE{2} - graddiv_E{2}/gam0^2 + curl_H{2}/gam0;
     KE{3} = KE{3} - graddiv_E{3}/gam0^2 + curl_H{3}/gam0; 
     KH{1} = KH{1} - graddiv_H{1}/gam0^2 - curl_E{1}/gam0;   
     KH{2} = KH{2} - graddiv_H{2}/gam0^2 - curl_E{2}/gam0;
     KH{3} = KH{3} - graddiv_H{3}/gam0^2 - curl_E{3}/gam0; 
end
function [KwE] = KopE(wE,input)
global nDIM;  

KwE = cell(1:nDIM); 
for n = 1:nDIM
    KwE{n} = Kop(wE{n},input.FFTG); 
end 
dummy = graddiv(KwE,input);           % dummy is temporary storage 
for n = 1:nDIM     
    KwE{n} = KwE{n} - dummy{n} / input.gamma_0^2; 
end        
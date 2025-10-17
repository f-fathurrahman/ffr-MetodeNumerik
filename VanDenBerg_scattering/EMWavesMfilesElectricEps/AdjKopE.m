function [Kf] = AdjKopE(f,input)
global nDIM;   
Kf = cell(1:nDIM);

for n = 1:nDIM
    Kf{n} = AdjKop(f{n},input.FFTG);
end

dummy = graddiv(Kf,input);            % dummy is temporary storage 

for n = 1:nDIM
    Kf{n} = Kf{n} - dummy{n} / conj(input.gamma_0^2);
end      
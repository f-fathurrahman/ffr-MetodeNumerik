function [Kv] = Kop(v,FFTG)
global nDIM;

Cv = zeros(size(FFTG));                            % make fft grid

if nDIM == 1; [N1,~]     = size(v);  Cv(1:N1,1)         = v; end  
if nDIM == 2; [N1,N2]    = size(v);  Cv(1:N1,1:N2)      = v; end  
if nDIM == 3; [N1,N2,N3] = size(v);  Cv(1:N1,1:N2,1:N3) = v; end 

Cv = fftn(Cv);   Cv = ifftn(FFTG .* Cv);           % convolution by fft

if nDIM == 1; Kv = Cv(1:N1,1);         end
if nDIM == 2; Kv = Cv(1:N1,1:N2);      end
if nDIM == 3; Kv = Cv(1:N1,1:N2,1:N3); end

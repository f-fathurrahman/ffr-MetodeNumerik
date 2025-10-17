function [Kf] = AdjKop(f,FFTG)
global nDIM;

Cf = zeros(size(FFTG));                            % make fft grid

if nDIM == 1; [N1,~]     = size(f);  Cf(1:N1,1)         = f; end  
if nDIM == 2; [N1,N2]    = size(f);  Cf(1:N1,1:N2)      = f; end  
if nDIM == 3; [N1,N2,N3] = size(f);  Cf(1:N1,1:N2,1:N3) = f; end  

Cf = fftn(Cf);   Cf = ifftn(conj(FFTG) .* Cf);     % convolution by fft

if nDIM == 1; Kf = Cf(1:N1,1);         end
if nDIM == 2; Kf = Cf(1:N1,1:N2);      end
if nDIM == 3; Kf = Cf(1:N1,1:N2,1:N3); end

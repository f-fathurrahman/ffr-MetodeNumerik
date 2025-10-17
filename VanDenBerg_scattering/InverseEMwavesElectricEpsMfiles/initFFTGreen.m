function input = initFFTGreen(input)
global nDIM;

if nDIM == 1                             % make one-dimensional FFT grid

   N1fft       = 2^ceil(log2(2*input.N1));
   x1(1:N1fft) = [0 : N1fft/2-1   N1fft/2 : -1 : 1] * input.dx; 
   temp.X1fft = x1';  % put it in single column  
 
elseif nDIM == 2                         % make two-dimensional FFT grid

   N1fft       = 2^ceil(log2(2*input.N1));
   N2fft       = 2^ceil(log2(2*input.N2));
   x1(1:N1fft) = [0 : N1fft/2-1   N1fft/2 : -1 : 1] * input.dx; 
   x2(1:N2fft) = [0 : N2fft/2-1   N2fft/2 : -1 : 1] * input.dx;
   [temp.X1fft,temp.X2fft] = ndgrid(x1,x2);  
 
elseif nDIM == 3                         % make three-dimensional FFT grid

   N1fft       = 2^ceil(log2(2*input.N1));
   N2fft       = 2^ceil(log2(2*input.N2));
   N3fft       = 2^ceil(log2(2*input.N3));
   x1(1:N1fft) = [0 : N1fft/2-1   N1fft/2 : -1 : 1] * input.dx; 
   x2(1:N2fft) = [0 : N2fft/2-1   N2fft/2 : -1 : 1] * input.dx;
   x3(1:N3fft) = [0 : N3fft/2-1   N3fft/2 : -1 : 1] * input.dx;
   [temp.X1fft,temp.X2fft,temp.X3fft] = ndgrid(x1,x2,x3);
   
end

% compute gam_0^2 * subdomain integrals  of Green function
  [IntG] = Green(temp,input);  

% apply n-dimensional Fast Fourier transform
  input.FFTG = fftn(IntG);


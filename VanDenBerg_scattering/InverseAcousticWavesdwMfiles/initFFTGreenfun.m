function [input] = initFFTGreenfun(input)

% make two-dimensional FFT grid
N1fft       = 2^ceil(log2(2*input.N1));
N2fft       = 2^ceil(log2(2*input.N2));
x1(1:N1fft) = [0 : N1fft/2-1   N1fft/2 : -1 : 1] * input.dx; 
x2(1:N2fft) = [0 : N2fft/2-1   N2fft/2 : -1 : 1] * input.dx;
[temp.X1fft,temp.X2fft] = ndgrid(x1,x2);
sign_x1 = [0 ones(1,N1fft/2-1)  0  -ones(1,N1fft/2-1)];
sign_x2 = [0 ones(1,N2fft/2-1)  0  -ones(1,N2fft/2-1)]; 
[Sign.X1fft,Sign.X2fft] = ndgrid(sign_x1,sign_x2);

% compute gam_0^2 * subdomain integrals of Green function
% and interface integrals  of  derivatives of Green function
  [IntG,IntdG]  = Greenfun(temp,Sign,input); 
  
% apply n-dimensional Fast Fourier transforms
  input.FFTG     = fftn(IntG);
  input.FFTdG{1} = fftn(IntdG{1});                                            
  input.FFTdG{2} = fftn(IntdG{2}); 
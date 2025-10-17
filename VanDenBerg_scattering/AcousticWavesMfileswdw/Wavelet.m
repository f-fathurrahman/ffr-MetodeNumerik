function [input] = Wavelet(input)

Nfft = 512;          fsamples = Nfft/2; 
df = 150 / fsamples;       dt = 1 / (df*Nfft);
  
% First derivative of Gaussian Wavelet using Matlab function gauswavf
  Nwavel = 21;
  [Waveltime,~] = gauswavf(-5,5,Nwavel,1); 

% Transform real function to frequency domain
  Wavelfrq = dt * conj(fft(Waveltime,Nfft));
  Wavelfrq(Nfft/2+1:Nfft) = 0;        % restrict to positive frequencies
  
% Transform to real function in time domain 
  Waveltime = 2*df * real(fft(Wavelfrq,Nfft)); 

figure; 
subplot(2,1,1); 
  plot((0:(fsamples-1))*df,abs(Wavelfrq(1:fsamples)),'b','LineWidth',1.5);  
  axis('tight');
  xlabel(' frequency  [Hz] \rightarrow ')
  ylabel(' |W(f)|  \rightarrow ');
  legend('frequency amplitude spectrum'); 
subplot(2,1,2); 
  plot((0:fsamples-1)*dt,Waveltime(1:fsamples),'r','LineWidth',1.5);
  axis('tight');
  xlabel(' time  [s] \rightarrow ') 
  ylabel(' W(t) \rightarrow ');
  legend('Wavelet in time domain')
  
input.Nfft      = Nfft;  
input.fsamples  = fsamples;
input.Wavelfrq = Wavelfrq(1:fsamples);
input.df        = df;
input.dt        = dt;

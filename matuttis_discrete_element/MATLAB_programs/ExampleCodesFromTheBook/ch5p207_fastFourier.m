% PURPOSE: Add/remove the last points in a 
%  periodic data sequence to experience which
%  data must be used to obtain a peak structure,
%  and see how the peaks broaden if too many
%  or not enough points are used
% CAVEAT: The last point 2*pi introduces noise,
%  it must be eliminated to get the 
%  "exact" Fourier spectrum.
% REVISION HISTORY: 21-May-2014 H.-G. Matuttis


clear all
format compact
l=50
x=linspace(0,2*pi,l);
%remove the end point (periodic):
% x=x(1:end-1);
% in this case the real part (middle window)
% will become numerically 0 (i.e. of the
% order of 10e-15, and of the imaginary part
% only the peaks (the wave number of the sine-curve)
% is left.
y=sin(x);
ffty=fft(y);
subplot(1,3,1)
plot(x,y,'*')
axis tight
subplot(1,3,2)
plot(real(ffty))
axis tight
subplot(1,3,3)
plot(imag(ffty))
axis tight
return
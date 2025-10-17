clear all; clc; close all; clear workspace
input = init();

%  (1) Compute analytically scattered field data --------------------------
       ForwardCanonicalObjects  

%  (2) Compute incident field ---------------------------------------------      
       u_inc = IncWave(input);
  
%  (3) Solve integral equation for wave field with CGFFT method -----------
tic;   u = ITERCGu(u_inc,input);
toc;   plotWavefield(u_inc,u,input); 
       w = input.CHI .* u;
       plotContrastSource(w,input);
      
%  (4) Compute scattered field data and plot fields and data --------------
       data = Dop(w,input);
       displayData(data,input);
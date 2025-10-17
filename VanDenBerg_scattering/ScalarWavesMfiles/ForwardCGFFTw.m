clear all; clc; close all; clear workspace
input = init();

%  (1) Compute analytically scattered field data --------------------------
       ForwardCanonicalObjects  

%  (2) Compute incident field ---------------------------------------------      
       u_inc = IncWave(input);
       
%  (3) Solve integral equation for contrast source with CGFFT method ------
       w =  ITERCGw(u_inc,input);             
       plotContrastSource(w,input);

%  (4) Compute synthetic data and plot fields and data --------------------
       data = Dop(w,input);
       displayData(data,input);    

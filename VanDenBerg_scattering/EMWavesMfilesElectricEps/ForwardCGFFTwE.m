clear all; clc; close all; clear workspace;
input = initEM();
       
%  (1) Compute analytically scattered field data --------------------------
       EMForwardCanonicalObjects
       
       plotEMcontrast(input); % plot permittivity / permeability contrast
       
%  (2) Compute incident field ---------------------------------------------      
       [E_inc,~] = IncEMwave(input);

%  (3) Solve integral equation for contrast source with CGFFT method ------
tic;   [w_E] = ITERCGwE(E_inc,input);
       plotContrastSourcewE(w_E,input);
toc
%  (4) Compute synthetic data and plot fields and data --------------------
       [Edata,Hdata] = DOPwE(w_E,input);             
       displayEdata(Edata,input);
       displayHdata(Hdata,input);
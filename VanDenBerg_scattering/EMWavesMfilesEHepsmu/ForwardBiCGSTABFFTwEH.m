clear all; clc; close all; clear workspace;
input = initEM();

%  (1) Compute analytically scattered field data --------------------------
       EMForwardCanonicalObjects
       
       plotEMcontrast(input); % plot permittivity / permeability contrast
       
%  (2) Compute incident field ---------------------------------------------      
       [E_inc,H_inc] = IncEMwave(input);        

% (3) Solve integral equation for contrast sources with FFT        
tic;  [w_E,w_H] = ITERBiCGSTABwEH(E_inc,H_inc,input);
toc      
% Plot electric and magnetic contrast sources
      plotContrastSourcewEH(w_E,w_H,input);
  
% (4) Compute synthetic data and plot fields and data ---------------------
      [Edata,Hdata] = DOPwEH(w_E,w_H,input);             
      displayEdata(Edata,input);
      displayHdata(Hdata,input);
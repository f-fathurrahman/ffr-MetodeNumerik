clear all; clc; close all; clear workspace;
input = initEM(); 
  
%  (1) Compute analytically scattered field data --------------------------
       EMForwardCanonicalObjects
       
       plotEMcontrast(input); % plot permittivity / permeability contrast
   
%  (2) Compute incident field ---------------------------------------------      
       [E_inc,~] = IncEMwave(input);

%  (3) Solve integral equation for contrast source with FFT ---------------
tic;   [w_E] = ITERBiCGSTABwE(E_inc,input);
       
      
toc;   plotContrastSourcewE(w_E,input);
      [E_sct] = KopE(w_E,input); 
       plotEtotalwavefield(E_inc,E_sct,input)

%  (4) Compute synthetic data and plot fields and data --------------------
       [Edata,Hdata] = DOPwE(w_E,input);             
       displayEdata(Edata,input);
       displayHdata(Hdata,input);
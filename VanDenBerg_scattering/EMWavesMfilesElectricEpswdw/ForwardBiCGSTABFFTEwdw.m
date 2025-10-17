clear all; clc; close all; clear workspace;
input = initEM();
       
%  (1) Compute analytically scattered field data --------------------------
       EMForwardCanonicalObjects
       
%      Input interface contrast 
       [input] = initEMcontrastIntf(input);

%  (2) Compute incident field ---------------------------------------------      
       [E_inc,dEinc] = IncEfield(input);   

%  (3) Solve integral equation for contrast source with FFT ---------------
tic;   [w_E,dwE] = ITERBiCGSTABEwdw(E_inc,dEinc,input);
toc        
%      Plot contrast and contrast sources
       plotContrastSourcewE(w_E,input); 
       plotEMInterfaceSourceE(dwE,input);

%  (4) Compute synthetic data and plot fields and data ---------------------
       [Edata,Hdata] = DOPEwdw(w_E,dwE,input);             
       displayEdata(Edata,input);
       displayHdata(Hdata,input);   
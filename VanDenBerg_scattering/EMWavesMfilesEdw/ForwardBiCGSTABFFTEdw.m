clear all; clc; close all; clear workspace;
input = initEM();

% (1) Compute alanalytically scattered field data
      EMForwardCanonicalObjects

%     Input interface contrast
      [input] = initEMcontrastIntf(input);

% (2) Compute incident field ----------------------------------------------      
      [Einc_tau] = IncEtaufield(input);  
      
% (3) Solve integral equation for contrast source with FFT ----------------
tic;  [dw] = ITERBiCGSTABEdw(Einc_tau,input);
toc
% Plot interface contrast and interface sources
      plotEMInterfaceSourceEdw(dw,input);

% (4) Compute synthetic data and plot fields and data ---------------------
      [Edata] = DOPmudw(dw,input);           
      displayEdata(Edata,input);
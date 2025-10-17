clear all; clc; close all; clear workspace
input = initAC();

%  (1) Comput analytically scattered field data
       AcForwardCanonicalObjects;

%      Input interface contrast 
       [input] = initAcousticContrastIntf(input);   

%  (2) Compute incident field at different grids      
       [p_inc,Pinc] = IncPressureWave(input);

%  (3) Solve contrast source integral equations with CGFFT method --------
tic;   [w,dw] = ITERCGwdw(p_inc,Pinc,input);
toc 
    %  Plot contrast sources, compute and plot the pressure wave fields
       plotContrastSource(w,input);
       plotInterfaceSource(dw,input);
     
% (4) Compute synthetic data and plot fields and data ---------------------
      data = DOPw(w,input);             
      data = DOPdw(dw,data,input);
      displayData(data,input);
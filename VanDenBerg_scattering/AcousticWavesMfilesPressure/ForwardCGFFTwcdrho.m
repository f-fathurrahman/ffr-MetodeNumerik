clear all; clc; close all; clear workspace
input = initAC();

%  (1) Compute analytically scattered field data --------------------------
       AcForwardCanonicalObjects

%      Input modified differentiable contrast functions -------------------    
       input = AcousticModifiedContrast(input); 

%  (2) Compute only the incident pressure field ---------------------------      
       [p_inc,~] = IncAcousticWave(input);  
       
%  (3) Solve contrast source integral equations with CGFFT method ---------
       input.method = 'analytical';   % choose 'analytical' or 'numerical'
       disp(['Derivatives are ' lower(input.method)]);
tic;   [w_crho,w_drho] = ITERCGwcdrho(p_inc,input); 
toc;    
    %  Plot contrast and contrast sources
       plotContrastSourceWcrho(w_crho,input);
       plotInterfaceSourceWdrho(w_drho,input);
 
%  (4) Compute synthetic data and plot fields and data --------------------
       data = DOPwcdrho(w_crho,w_drho,input);             
       displayData(data,input);
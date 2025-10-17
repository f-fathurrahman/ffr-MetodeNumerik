clear all; clc; close all; clear workspace
input = initAC();

%  (1) Compute analytically scattered field data --------------------------
       AcForwardCanonicalObjects
       
%  (2) Compute incident field ---------------------------------------------      
       [p_inc,Zv_inc] = IncAcousticWave(input);

%  (3) Solve contrast source integral equations with FFT method -----------
       input.method = 'analytical';   % choose 'analytical' or 'numerical'
       disp(['Derivatives are ' lower(input.method)]); 
  tic;     [w_p,w_Zv] = ITERBiCGSTABwpZv(p_inc,Zv_inc,input);      toc;
         
    %  Plot contrast sources, compute and plot the pressure wave fields
       plotContrastSourcewp(w_p,input);
       plotContrastSourcewZv(w_Zv,input);
       [Kf,~] = KOPpZv(w_p,w_Zv,input);  p = p_inc + Kf;
       plotPressureWavefield(p_inc,p,input);
     
% (4)  Compute synthetic data and plot fields and data --------------------
       data = DOPwpZv(w_p,w_Zv,input);             
       displayData(data,input);           
clear all; clc; close all; clear workspace
input = initAC();

% (1) Compute analytically scattered field data --------------------------
      disp('Running AcousticDataSctCircle');  AcousticDataSctCircle;
      
      %load DataAnalytic; 
      %load DataIntEq;   dataCircle = data;       

data = zeros(input.NR,input.NS);
for q = 1 : input.NS   
% (2) Compute incident field at different grids      
      Pinc = IncPressureWave(q,input);

% (3) Solve contrast source integral equations with CGFFT method ---------
      % dw = ITERCGdw(Pinc,input); 
        dw = BICGSTABdw(Pinc,input);
   
% (4) Compute synthetic data and plot fields and data ---------------------
      data(:,q) = DOPMdw(dw,input);     
end % q_loop 
  save DataIntEq data;
            
 % Plot data at a number of receivers -------------------
   figure(10);
   subplot(1,2,1) 
     PlotData(dataCircle,input);
     title('|p^{sct}(x_p^R,x_q^S)|=|data_{Bessel}|')
   subplot(1,2,2)
     PlotData(data-dataCircle,input);
     
 % Compute normalized difference of norms  
   error = num2str(norm(data(:)-dataCircle(:),1)/norm(dataCircle(:),1));
   disp(['error=' error]);      title('|data_{Int.Eq.} - data_{Bessel}|')
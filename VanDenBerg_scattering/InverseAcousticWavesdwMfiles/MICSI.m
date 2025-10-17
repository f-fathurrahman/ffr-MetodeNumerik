clear all; clc; close all; clear workspace
input = initAC();  

% (1) Compute Exact data using analytic Bessel function expansion --------
      disp('Running AcousticDataSctCircle');  AcousticDataSctCircle;
      
      load DataAnalytic; 
      % load DataIntEq;   dataCircle = data;         % Inverse Crime
      
      % add complex-valued noise to data  
      if input.Noise == 1
         dataCircle = dataCircle + mean(abs(dataCircle(:))) * input.Rand;
      end
         
%  (2) Apply MRCCSI method -----------------------------------------------
       [dw] = ITERCGMDdw(dataCircle,input); 
       
%  (3) Compute model data explained by reconstructed contrast sources ----
       data = zeros(input.NR,input.NS);
       for q = 1 : input.NS
         dwq{1} = dw{1,q};
         dwq{2} = dw{2,q};
         data(:,q) = DOPMdw(dwq,input); 
       end % q_loop 

%   plot data and compare to synthetic data
    figure(20);
    subplot(1,2,1) 
      PlotData(dataCircle,input);
      title('|p^{sct}(x_p^R,x_q^S)|=|data_{Bessel}|')
    subplot(1,2,2)
      PlotData(data-dataCircle,input);
%   Compute normalized difference of norms  
    error = num2str(norm(data(:)-dataCircle(:),1)/norm(dataCircle(:),1));
    disp(['error=' error]);      title('|data_{Int.Eq.} - data_{Exact}|')
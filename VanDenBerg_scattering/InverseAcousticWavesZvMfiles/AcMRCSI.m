clear all; clc; close all; clear workspace
input = initAC();

% (1) Compute Exact data using analytic Bessel function expansion --------
      disp('Running AcousticDataSctCircle');   AcousticDataSctCircle;
      
      load DataAnalytic;
      
      % load DataIntEq;   dataCircle = data;        % Inverse Crime

      % add complex-valued noise to data  
      if input.Noise == 1
         dataCircle = dataCircle + mean(abs(dataCircle(:))) * input.Rand;
      end      
      
%  (2) Compute Green functions for sources and receivers ------------------
       [~,dG_S,~,dG_R] = dGreenSourcesReceivers(input);
        
%  (3) Apply MRCSI method -------------------------------------------------
       Zvinc = cell(2,input.NS);
       for q = 1 : input.NS  
          Zvinc{1,q} = - (1/input.gamma_0) * dG_S{1,q};
          Zvinc{2,q} = - (1/input.gamma_0) * dG_S{2,q};
       end % q_loop 
       clear G_S dG_S 
       
       [wZv] = ITERCGMDwZv(Zvinc,dG_R,dataCircle,input);

%  (4) Compute model data explained by reconstructed contrast souces ------
       wZvq = cell(2,input.NS);
        data = zeros(input.NR,input.NS);
       for q = 1 : input.NS
           wZvq{1} = wZv{1,q};
           wZvq{2} = wZv{2,q};
         data(:,q) = DOPMwZv(dG_R,wZvq,input); 
       end % q_loop 

%  Plot data and compare to synthetic data
   figure(9);
   subplot(1,2,1);   
      PlotData(dataCircle,input);
      title('|p^{sct}(x_p^R,x_q^S)|=|data_{Bessel}|')
   subplot(1,2,2);   
      PlotData(data-dataCircle,input);
      title('|data_{Int.Eq.} - data_{Exact}|')
      
%  Compute normalized difference of norms  
   error = num2str(norm(data(:)-dataCircle(:),1)/norm(dataCircle(:),1));
   disp(['error=' error]);      
clear all; clc; close all; clear workspace
input = initAC();

%  (1) Compute analytically scattered field data --------------------------
       disp('Running AcousticDataSctCircle');  AcousticDataSctCircle;
             
%  (2) Compute Green functions for sources and receivers ------------------
       [G_S,dG_S,G_R,dG_R] = dGreenSourcesReceivers(input);
       
%  (3) Solve contrast source integral equations and compute sunthetic data
       data = zeros(input.NR,input.NS);
       w_Zv = cell(2,input.NS);
       for q = 1 : input.NS  
          Zv_inc{1} = - (1/input.gamma_0) * dG_S{1,q};
          Zv_inc{2} = - (1/input.gamma_0) * dG_S{2,q};
         %[wZv]     = ITERCGwZv(Zv_inc,input);
          [wZv]     = BICGSTABwZv(Zv_inc,input);
          data(:,q) = DOPMwZv(dG_R,wZv,input);
       end % q_loop 
       
       save DataIntEq data;
                 
%  (4) Plot data at a number of receivers -------------------
       figure(5);
       subplot(1,2,1) 
          PlotData(dataCircle,input);
          title('|p^{sct}(x_p^R,x_q^S)|=|data_{Bessel}|')
       subplot(1,2,2)
          PlotData(data-dataCircle,input);
          title('|data_{Int.Eq.} - data_{Bessel}|');
     % Compute normalized difference of norms  
       error=num2str(norm(data(:)-dataCircle(:),1)/norm(dataCircle(:),1));
       disp(['error=' error]);    
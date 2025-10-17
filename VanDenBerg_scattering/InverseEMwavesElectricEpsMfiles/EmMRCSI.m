clear all; clc; close all; clear workspace;
input = initEM();
%  (1) Compute Exact data using anlaytic Bessel function expansion -------
       disp('Running EMdataCircle'); EMdataCircle;
       load Edata E1data E2data;
     % load EMdataIntEq  E1DATAnum  E2DATAnum;      
     % E1data = E1DATAnum;  E2data = E2DATAnum;          % Inverse Crime
     
     % add complex-valued noise to data 
       if input.Noise == 1
          E1data = E1data + mean(abs(E1data(:))) * input.Rand;
          E2data = E2data + mean(abs(E2data(:))) * input.Rand;
       end      
     
%  (2) Compute Green function for sources and receivers -------------------     
       pG = EGreenSourceReceivers(input);
       
%  (3) Apply MRCCSI method and plot inverted permittivity contrast --------
       Einc = cell(2,input.NS);
       for q = 1 : input.NS  
          [E_incq]   = IncEMwave(pG,q,input);
           Einc{1,q} = E_incq{1};   Einc{2,q} = E_incq{2};
       end %q_loop
       [wE] = ITERMDwE(Einc,pG,E1data,E2data,input);      
           
%  (4) Compute model data explained by reconstructed contrast sources -----
       E1dataInv= zeros(input.NR,input.NS);  
       E2dataInv= zeros(input.NR,input.NS);
       for q = 1 : input.NS
         w_E{1} = wE{1,q};  w_E{2} = wE{2,q};
         [Ercv] = DOPwE(pG,w_E,input); 
          E1dataInv(:,q) = Ercv{1};   E2dataInv(:,q) = Ercv{2};
       end % q_loop  
       
    % Plot data and compare to synthetic data
      figure(9);                               
          subplot(1,2,1);  PlotData(E1data,input);                     
          title('|E_1^{sct}(x_p^R,x_q^S)| = |E_1data_{Bessel}|');
          subplot(1,2,2);  PlotData(E1dataInv-E1data,input); 
          title('|E_1data_{Int.Eq.} - E_1data_{Bessel}|'); 
      figure(10);                             
          subplot(1,2,1);  PlotData(E2data,input);
          title('|E_2^{sct}(x_p^R,x_q^S)| = |E_2data_{Bessel}|');
          subplot(1,2,2);  PlotData(E2dataInv-E2data,input);             
          title('|E_2data_{Int.Eq.} - E_2data_{Bessel}|');       
   % Compute normalized difference of norms   
     error1 = num2str(norm(E1dataInv(:)-E1data(:),1)/norm(E1data(:),1));  
     error2 = num2str(norm(E2dataInv(:)-E2data(:),1)/norm(E2data(:),1));
     disp(['error1 =' error1]);  disp(['error2 =' error2]);
     
     save MRCSIinv E1dataInv E2dataInv; 
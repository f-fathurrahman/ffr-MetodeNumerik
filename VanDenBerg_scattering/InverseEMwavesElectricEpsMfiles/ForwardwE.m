clear all; clc; close all; clear workspace;
input = initEM();
       
%  (1) Compute analytically scattered field data --------------------------
       disp('Running EMdataCircle'); EMdataCircle;
       
       plotEMcontrast(input); % plot permittivity / permeability contrast
       
%  (2) Compute Green function for sources and receivers -------------------     
       pG = EGreenSourceReceivers(input);
    
%  (3) Solve integral equation for contrast source with CGFFT method ------        
       E1dataNum = zeros(input.NR,input.NS);  
       E2dataNum = zeros(input.NR,input.NS);
       for q = 1 : input.NS  
           [E_incq] = IncEMwave(pG,q,input); 
            %[w_E]  = ITERCGwE(E_incq,input);
             [w_E]  = BICGSTABwE(E_incq,input);
           [Ercv]   = DOPwE(pG,w_E,input); 
           E1dataNum(:,q) = Ercv{1}; E2dataNum(:,q) = Ercv{2};      
       end %q_loop
       
       save EMdataIntEq  E1dataNum  E2dataNum;    
       
%  (4) Plot data at receiver locations ----------------------------------    
       figure(9);                              
          subplot(1,2,1);  PlotData(E1data,input);                     
          title('|E_1^{sct}(x_p^R,x_q^S)| = |E_1data_{Bessel}|');
          subplot(1,2,2);  PlotData(E1dataNum-E1data,input); 
          title('|E_1data_{Int.Eq.} - E_1data_{Bessel}|'); 
       figure(10);                            
          subplot(1,2,1);  PlotData(E2data,input);
          title('|E_2^{sct}(x_p^R,x_q^S)| = |E_2data_{Bessel}|');
          subplot(1,2,2);  PlotData(E2dataNum-E2data,input);             
          title('|E_2data_{Int.Eq.} - E_2data_{Bessel}|');  
    
     % Compute normalized difference of norms      
       error1 = num2str(norm(E1dataNum(:)-E1data(:),1)/norm(E1data(:),1));  
       disp(['error1 =' error1]); 
       error2 = num2str(norm(E2dataNum(:)-E2data(:),1)/norm(E2data(:),1)); 
       disp(['error2 =' error2]);  
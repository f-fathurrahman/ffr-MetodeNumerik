clear all; clc; close all; clear workspace;
input = initEM();

%  (1) Compute analytically scattered field data for monitor model -------
       disp('Running EMdataCircle'); EMdataCircle;    
       
%      Input interface contrast 
       [input] = initEMcontrastIntf(input); 
          
%  (2) Compute Green function for sources and receivers ------------------     
       pG = EGreenSourceReceivers(input);

%  (3) Compute base-model data using either BICGSTABwE or BICGSTABwdwE ---
       E1DATAwE = zeros(input.NR,input.NS);    
       E2DATAwE = zeros(input.NR,input.NS);
       for q = 1 : input.NS
           [E_incq,~]      = IncEMdwave(pG,q,input); 
               [w_E]       = BICGSTABwE(E_incq,input); 
               [ErSct]     = DOPwE(pG,w_E,input);         
           E1DATAwE(:,q)   = ErSct{1};
           E2DATAwE(:,q)   = ErSct{2}; 
       end %q_loop
       
       E1DATAwdwE = zeros(input.NR,input.NS);    
       E2DATAwdwE = zeros(input.NR,input.NS);
       for q = 1 : input.NS
           [E_incq,dEincq] = IncEMdwave(pG,q,input);
               [w_E,dwE]   = BICGSTABwdwE(E_incq,dEincq,input);
               [ErSct]     = DOPEwdw(w_E,dwE,input);      
           E1DATAwdwE(:,q) = ErSct{1};
           E2DATAwdwE(:,q) = ErSct{2};    
       end %q_loop

%  (4) Plot monitor-model and base-model data using wE or wdwE method ----
       [difE1,difE2] =  ...
           PlotDataDif(E1data,E2data,E1DATAwE,E2DATAwE,input);       
       disp(['difwE1 =' difE1]); 
       disp(['difwE2 =' difE2]);  
       
       [difE1,difE2] =  ...
           PlotDataDif(E1data,E2data,E1DATAwdwE,E2DATAwdwE,input);       
       disp(['difwdwE1 =' difE1]); 
       disp(['difwdwE2 =' difE2]);
       
%  (5) Focus the two type of data sets using wE or wdwE method ---------- 
       Focusing(E1data,E2data,E1DATAwE,  E2DATAwE,  input);    
       Focusing(E1data,E2data,E1DATAwdwE,E2DATAwdwE,input); 
       
       PlotDefinitionModels
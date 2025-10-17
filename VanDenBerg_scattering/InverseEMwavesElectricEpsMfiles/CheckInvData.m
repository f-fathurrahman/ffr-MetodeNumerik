clear all; clc; close all; clear workspace;
input = initEM();
%     Input interface contrast 
      [input] = initEMcontrastIntf(input); 
%     
%     Load the analytic data and the inverted data from MRCSI ------------
      load Edata E1data E2data;
      load MRCSIinv E1dataInv E2dataInv; 
      load Contrast256 CHI_256; 
      input.CHI_eps = CHI_256;  
      x1 = input.X1(:,1);   x2 = input.X2(1,:);
      figure(5); IMAGESC(x1,x2,real(input.CHI_eps)); 
      title('\fontsize{13} Reconstructed Contrast Re[\chi^\epsilon] ');       

% (1) Compute normalized difference of norms   
      error1 = num2str(norm(E1dataInv(:)-E1data(:),1)/norm(E1data(:),1));  
      error2 = num2str(norm(E2dataInv(:)-E2data(:),1)/norm(E2data(:),1));
      disp(['inv ErrorE1 =' error1]);
      disp(['inv ErrorE2 =' error2]);
   
%     Compute Green function for sources and receivers -------------------     
      pG = EGreenSourceReceivers(input);   
%     Solve integral equation using wE-method and wdwE-method ------------    
      E1DATAwE = zeros(input.NR,input.NS);    
      E2DATAwE = zeros(input.NR,input.NS);    
      E1DATAwdwE = zeros(input.NR,input.NS); 
      E2DATAwdwE = zeros(input.NR,input.NS);    
      for q = 1 : input.NS
           [E_incq,dEincq] = IncEMdwave(pG,q,input); 
               [w_E]       = BICGSTABwE(E_incq,input); 
               [ErSct]     = DOPwE(pG,w_E,input);         
           E1DATAwE(:,q)   = ErSct{1};
           E2DATAwE(:,q)   = ErSct{2};
               [w_E,dwE]   = BICGSTABwdwE(E_incq,dEincq,input);
               [ErSct]     = DOPEwdw(w_E,dwE,input);      
           E1DATAwdwE(:,q) = ErSct{1};
           E2DATAwdwE(:,q) = ErSct{2};    
      end %q_loop  

% (2) Compute normalized difference of norms (wE-method)  
      error1 = num2str(norm(E1DATAwE(:)-E1data(:),1)/norm(E1data(:),1));  
      error2 = num2str(norm(E2DATAwE(:)-E2data(:),1)/norm(E2data(:),1));
      disp(['wE ErrorE1 =' error1]);  
      disp(['wE ErrorE2 =' error2]);
     
% (3) Compute normalized difference of norms (wdE-method)  
      error1 = num2str(norm(E1DATAwdwE(:)-E1data(:),1)/norm(E1data(:),1));  
      error2 = num2str(norm(E2DATAwdwE(:)-E2data(:),1)/norm(E2data(:),1));
      disp(['wdwE Error1 =' error1]);  
      disp(['wdwE error2 =' error2]);
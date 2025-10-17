clear all; clc; close all; clear workspace
input = init();

%  (1) Compute analytically scattered field data --------------------------
       disp('Running DataSctCircle');  DataSctCircle;
       
%  (2) Compute Green function for sources and receivers -------------------     
       [G_S,G_R] = GreenSourcesReceivers(input);

%  (3) Solve contrast source integral equation and compute synthetic data        
       data = zeros(input.NR,input.NS);
          w = cell(1,input.NS);            
       for q = 1 : input.NS  
   %       w{q} = ITERCGw(G_S{q},input);   % Note u_inc(:,:) = G_S{q}(:,:)            
           w{q} = BICGSTABw(G_S{q},input);   
         data(:,q) = DopM(G_R,w{q},input); 
       end % q_loop 
       save DataIntEq data;

%  (4) Plot data at a number of receivers ---------------------------------
       figure(3);
       subplot(1,2,1) 
          PlotData(dataCircle,input);
          title('|u^{sct}(x_p^R,x_q^S)|=|data_{Bessel}|')
       subplot(1,2,2)
          PlotData(data-dataCircle,input);  
          title('|data_{Int.Eq.} - data_{Bessel}|') 
  
     % Compute normalized difference of norms  
     error = num2str(norm(data(:)-dataCircle(:),1)/norm(dataCircle(:),1));
     disp(['error=' error]);      
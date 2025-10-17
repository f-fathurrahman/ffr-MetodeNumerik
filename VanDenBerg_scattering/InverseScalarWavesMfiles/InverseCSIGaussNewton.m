clear all; clc; close all; clear workspace
input = init();

% (1) Compute Exact data using analytic Bessel function expansion -------
       disp('Running DataSctCircle');  DataSctCircle;  
       load DataAnalytic; 
       %load DataIntEq;   dataCircle = data;
       % add complex-valued noise to data  
       if input.Noise == 1
          dataCircle = dataCircle + mean(abs(dataCircle(:))) * input.Rand;
       end      
           
% (2) Compute Green function for sources and receivers ------------------
      [G_S,G_R] = GreenSourcesReceivers(input);  
       
% (3) Apply MRCCSI method -----------------------------------------------
      [w,CHI] = ITERCGMDGaussNewton(G_S,G_R,dataCircle,input);  
      
% Figure of reconstructed contrast (analytic form)
  figure(5);  IMAGESC(input.X1(:,1),input.X2(1,:),real(CHI)); 
  figure(6);  mesh(abs(CHI),'edgecolor', 'k'); 
              view(37.5,45); axis xy; axis('off'); axis('tight') 
              axis([1 input.N1 1 input.N2  0 1.5*max(abs(input.CHI(:)))]);
              
% Figure of reconstructed contrast parameters
  DisplayParameters;     
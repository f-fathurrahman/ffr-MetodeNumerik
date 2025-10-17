clear all; clc; close all; clear workspace
input = init();

% (1) Compute Exact data using analytic Bessel function expansion --------
      disp('Running DataSctCircle');  DataSctCircle;  
      
      load DataAnalytic; 
      % load DataIntEq;   dataCircle = data;       % Inverse Crime
       
      % add complex-valued noise to data  
      if input.Noise == 1
         dataCircle = dataCircle + mean(abs(dataCircle(:))) * input.Rand;
      end     
       
% (2) Compute Green function for sources and receivers -------------------
      [G_S,G_R] = GreenSourcesReceivers(input);  

% (3) Apply CSI method ---------------------------------------------------
      [w,CHI] = ITERCGMD(G_S,G_R,dataCircle,input);  
      
   % Figure of contrast as image and surface plot 
     x1 = input.X1(:,1);   x2 = input.X2(1,:);
     figure(5); IMAGESC(x1,x2,real(CHI)); 
  
     figure(6); mesh(abs(CHI),'edgecolor', 'k'); 
                view(37.5,45); axis xy; axis('off'); axis('tight') 
                axis([1 input.N1 1 input.N2 0 1.5*max(abs(input.CHI(:)))])
   
% (4) Compute model data explained by reconstructed contrast souces ------
      data = zeros(input.NR,input.NS);
      for q = 1 : input.NS
         data(:,q) = DopM(G_R,w{q},input);
      end % q_loop 
       
   % plot data and compare to synthetic data
     error = num2str(norm(data(:)-dataCircle(:),1)/norm(dataCircle(:),1));
     disp(['data error=' error]);  
     figure(7);
     subplot(1,2,1) 
       imagesc(1:input.NR,1:input.NS, abs(data)); 
       xlabel('x_p^R \rightarrow'); 
       ylabel('x_q^S \rightarrow');  
       axis('equal','tight'); axis xy; colorbar('hor'); colormap jet;
       title('|data_{Exact}|')
     subplot(1,2,2)
       imagesc(1:input.NR,1:input.NS, abs(data-dataCircle));
       xlabel('x_p^R \rightarrow'); 
       ylabel('x_q^S \rightarrow');  
       axis('equal','tight'); axis xy; colorbar('hor'); colormap jet;
       title('|data_{Inverted} - data_{Exact}|')  
function Focusing(E1monitor,E2monitor,E1base,E2base,input)
x1 = input.X1(:,1);   x2 = input.X2(1,:);   

% Plot base model %-------------------------------------------------------
figure; 
subplot(1,3,1);  
  IMAGESC(x1,x2,real(input.CHI_eps));                 
  grid on; ax = gca; ax.GridColor = 'w'; ax.GridAlpha =0.5;
  title('\fontsize{11}Base model');  xlabel(''); ylabel(''); 
 
% Focus E1 and E2 data with phase information %---------------------------
  [dataFocus] = divFocus(E1monitor-E1base,E2monitor-E2base,input);
  
subplot(1,3,2);
   IMAGESC(x1,x2, abs(dataFocus)); 
   grid on; ax = gca; ax.GridColor = 'w';  ax.GridAlpha =0.5;
   title('\fontsize{11}Data with phase');
   xlabel(''); ylabel(''); 
   hold on; phi = 0:.01:2*pi; 
   plot(input.a*sin(phi),input.a*cos(phi),'-y','LineWidth',1.5); hold off;
   colormap(jet.^(1/4));
   h = caxis;  margin = (h(2)-h(1))/4;
   caxis([h(1)+ margin; h(2)-margin]);
       
% Focus Amplitude E1 and E2 data %----------------------------------------
  fase1base = E1base./abs(E1base); 
  fase2base = E2base./abs(E2base);
  % add fase factors to amplitudes of datawE
    E1monitor = abs(E1monitor) .* fase1base;
    E2monitor = abs(E2monitor) .* fase2base;
    [dataFocus] = divFocus(E1monitor-E1base,E2monitor-E2base,input);   
    
subplot(1,3,3)
   IMAGESC(x1,x2, abs(dataFocus));  
   title('\fontsize{11}Amplitude data');
   grid on; ax = gca; ax.GridColor = 'w';  ax.GridAlpha =0.5;
   xlabel(''); ylabel(''); 
   hold on; phi = 0:.01:2*pi; 
   plot(input.a*sin(phi),input.a*cos(phi),'-y','LineWidth',1.5); hold off; 
   colormap(jet.^(1/4));
   h = caxis;  margin = (h(2)-h(1))/4;
   caxis([h(1)+ margin; h(2)-margin]);
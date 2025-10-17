function  [difE1,difE2] =  ...
              PlotDataDif(E1monitor,E2monitor,E1base,E2base,input)

absEmonitor = sqrt(abs(E1monitor).^2 + abs(E2monitor).^2);
absEbase    = sqrt(abs(E1base).^2 + abs(E2base).^2); 

set(figure,'Units','centimeters','Position', [0 0 25 25] );
    sp = subplot(1,3,1);
         IMAGESC(input.X1(:,1),input.X1(:,1),(input.CHI_eps));                 
         grid on; ax = gca; ax.GridColor = 'w'; ax.GridAlpha =0.5;
         title('\fontsize{12}Base model'); 
         pos = get(sp,'Position'); 
         pos = pos+[0.05 0.05 -0.05 -0.057]; set(sp,'Position',pos);
    subplot(1,3,2);    
         PlotData(absEbase,input);                    
         title('\fontsize{12}|Edata|_{base}');        
    subplot(1,3,3);  
         PlotData(absEmonitor-absEbase,input);
         title('\fontsize{12}|Edata|_{monitor}-|Edata|_{base}');      
         colormap(jet.^(1/4));
         h = caxis;  margin = (h(2)-h(1))/4; 
         caxis([h(1)+margin; h(2)-margin]);        

% Compute normalized difference of norms      
     difE1 = num2str(norm(E1monitor(:)-E1base(:),1)/norm(E1monitor(:),1));  
     difE2 = num2str(norm(E2monitor(:)-E2base(:),1)/norm(E2monitor(:),1)); 
 
function  displayEdata(Edata,input)
global nDIM;

if nDIM == 2;  % plot data at a number of receivers -----------------------
    
   if exist(fullfile(cd,'EDATA2D.mat'), 'file');  
      load EDATA2D Edata2D; 
      E_error = num2str(norm(Edata(:)-Edata2D(:),1)/norm(Edata2D(:),1));    
      disp(['E-error=' E_error]);
   end
   set(figure,'Units','centimeters','Position', [5 5 18 7]);
   angle = input.rcvr_phi * 180 / pi;
   
   if exist(fullfile(cd,'EDATA2D.mat'), 'file')
      plot(angle,abs(Edata),'--r',angle,abs(Edata2D),'b')
      legend('Integral-equation method', ...
              'Bessel-function method','Location','Best');
      text(50,0.8*max(abs(Edata)), ...
      ['Error(E^{sct}) = ' E_error '  '],'EdgeColor','red','Fontsize',11);
   else
      plot(angle,abs(Edata),'b')
      legend('Bessel-function method','Location','Best');
   end 
   title('\fontsize{12} scattered E data in 2D');   axis tight;
   xlabel('observation angle in degrees'); 
   ylabel('abs(data) \rightarrow');  
    
elseif nDIM == 3;  % plot data at a number of receivers ------------------- 
    
   if exist(fullfile(cd,'EDATA3D.mat'), 'file');  
      load EDATA3D Edata3D ;  
      E_error = num2str(norm(Edata(:)-Edata3D(:),1)/norm(Edata3D(:),1));
      disp(['E-error=' E_error]);    
   end
   set(figure,'Units','centimeters','Position', [5 5 18 7]);
   angle = input.rcvr_phi * 180 / pi;
   
   if exist(fullfile(cd,'EDATA3D.mat'), 'file')
      plot(angle,abs(Edata),'--r',angle,abs(Edata3D),'b')
      legend('Integral-equation method', ...
              'Bessel-function method','Location','Best');
      text(50,0.8*max(abs(Edata)), ...
      ['Error(E^{sct}) = ' E_error '  '],'EdgeColor','red','Fontsize',11);   
   else
      plot(angle,abs(Edata),'b')
      legend('Bessel-function method','Location','Best');
   end 
   title('\fontsize{12} scattered E data in 3D');  axis tight;
   xlabel('observation angle in degrees'); 
   ylabel('abs(data) \rightarrow');    
end % if
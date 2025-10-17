function  displayHdata(Hdata,input)
global nDIM;

if nDIM == 2;  % plot data at a number of receivers -----------------------
    
   if exist(fullfile(cd,'HDATA2D.mat'), 'file');  
     load HDATA2D Hdata2D; 
      H_error = num2str(norm(Hdata(:)-Hdata2D(:),1)/norm(Hdata2D(:),1));    
      disp(['H-error=' H_error]);
   end
   set(figure,'Units','centimeters','Position', [5 5 18 7]);
   angle = input.rcvr_phi * 180 / pi;
   
   if exist(fullfile(cd,'HDATA2D.mat'), 'file')
      plot(angle,abs(Hdata),'--r',angle,abs(Hdata2D),'b')
      legend('Integral-equation method', ...
              'Bessel-function method','Location','Best');
      text(50,0.8*max(abs(Hdata)), ...
      ['Error(Z_0H^{sct}) = ' H_error '  '],'EdgeColor','red','Fontsize',11);
   else
      plot(angle,abs(Hdata),'b')
      legend('Bessel-function method','Location','Best');
   end 
   title('\fontsize{12} scattered Z_0H data in 2D');   axis tight;
   xlabel('observation angle in degrees'); 
   ylabel('abs(data) \rightarrow');    
     
elseif nDIM == 3;  % plot data at a number of receivers -------------------
    
   if exist(fullfile(cd,'HDATA3D.mat'), 'file');  
     load HDATA3D Hdata3D; 
     H_error = num2str(norm(Hdata(:)-Hdata3D(:),1)/norm(Hdata3D(:),1));
     disp(['H-error=' H_error]);    
   end
   set(figure,'Units','centimeters','Position', [5 5 18 7]);
   angle = input.rcvr_phi * 180 / pi;
   
   if exist(fullfile(cd,'HDATA3D.mat'), 'file')
      plot(angle,abs(Hdata),'--r',angle,abs(Hdata3D),'b')
      legend('Integral-equation method', ...
              'Bessel-function method','Location','Best');
      text(50,0.8*max(abs(Hdata)), ...
      ['Error(Z_0H^{sct}) = ' H_error '  '],'EdgeColor','red','Fontsize',11);        
   else
      plot(angle,abs(Hdata),'b')
      legend('Bessel-function method','Location','Best');
   end 
   title('\fontsize{12} scattered Z_0H data in 3D');   axis tight;
   xlabel('observation angle in degrees'); 
   ylabel('abs(data) \rightarrow');    
end % if
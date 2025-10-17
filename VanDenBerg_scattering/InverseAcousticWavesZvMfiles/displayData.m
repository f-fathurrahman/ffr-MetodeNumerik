function  displayData(data,input)
global nDIM;

if nDIM == 1     % display scattered data at two receivers
  disp(['scattered wave amplitude above slab = ' num2str(abs(data(1)))]);
  disp(['scattered wave amplitude below slab = ' num2str(abs(data(2)))]); 
  if exist(fullfile(cd, 'data1D.mat'), 'file');     load data1D data1D;   
     disp(['analytic data above slab = ' num2str(abs(data1D(1)))]);
     disp(['analytic data below slab = ' num2str(abs(data1D(2)))]);
     disp(['error=' num2str(norm(data(:)-data1D(:),1)/norm(data1D(:),1))]);
  end

elseif nDIM == 2  % plot data at a number of receivers -------------------
  if exist(fullfile(cd,'DATA2D.mat'), 'file');      load DATA2D data2D; 
     error = num2str(norm(data(:)-data2D(:),1)/norm(data2D(:),1));
     disp(['error=' error]);  
  end   
  set(figure,'Units','centimeters','Position', [5 5 18 7]);
  angle = input.rcvr_phi * 180 / pi;
  if exist(fullfile(cd,'DATA2D.mat'), 'file')
       plot(angle,abs(data),'--r',angle,abs(data2D),'b')
       legend('Integral-equation method', ...
              'Bessel-function method','Location','NorthEast');
       text(50,0.8*max(abs(data)), ...
            ['Error^{sct} = ' error '  '],'EdgeColor','red','Fontsize',11);
  else plot(angle,abs(data),'b')
       legend('Bessel-function method','Location','NorthEast');
  end 
  title('\fontsize{12} scattered wave data in 2D');   axis tight;
  xlabel('observation angle in degrees'); ylabel('abs(data) \rightarrow');  
  
elseif nDIM == 3  % plot data at a number of receivers ------------------- 
  if exist(fullfile(cd,'DATA3D.mat'), 'file');      load DATA3D data3D;  
     error = num2str(norm(data(:)-data3D(:),1)/norm(data3D(:),1));   
     disp(['error=' error]);  
  end
  set(figure,'Units','centimeters','Position', [5 5 18 7]);
  angle = input.rcvr_phi * 180 / pi;
  if exist(fullfile(cd,'DATA3D.mat'), 'file')
       plot(angle,abs(data),'--r',angle,abs(data3D),'b')
       legend('Integral-equation method', ...
              'Bessel-function method','Location','Best');
       text(50,0.8*max(abs(data)), ...
            ['Error^{sct} = ' error '  '],'EdgeColor','red','Fontsize',11);
  else plot(angle,abs(data),'b')
       legend('Bessel-function method','Location','Best');
  end   
  title('\fontsize{12} scattered wave data in 3D');   axis tight;
  xlabel('observation angle in degrees'); ylabel('abs(data) \rightarrow'); 
end % if
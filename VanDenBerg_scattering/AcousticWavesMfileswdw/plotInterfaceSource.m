  function plotInterfaceSource(dw,input)
global nDIM;  dx =input.dx; Rfl = input.Rfl;

if nDIM == 1           % Plot 1D contrast/source distribution -------------
    x1 = input.X1; 
    set(figure,'Units','centimeters','Position',[5 5 18 12]);
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,-30:30,Rfl{1});  
       title('\fontsize{13} R^{(1)}=(\rho^+-\rho^-) / (\rho^++\rho^-)');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,-30:30, abs(dw{1}));   
       title('\fontsize{13} abs(dw^{(1)})');  
 
elseif  nDIM == 2     % Plot 2D contrast/source distribution -------------
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2,Rfl{1});
       title('\fontsize{13} R^{(1)}=(\rho^+-\rho^-) / (\rho^++\rho^-)');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,x2, abs(dw{1}));   
       title('\fontsize{13} abs(dw^{(1)})');
    set(figure,'Units','centimeters','Position',[5 5 18 12]);       
    subplot(1,2,1);  
       IMAGESC(x1,x2+dx/2,Rfl{2});
       title('\fontsize{13} R^{(2)}=(\rho^+-\rho^-) / (\rho^++\rho^-)');
    subplot(1,2,2);  
       IMAGESC(x1,x2+dx/2,abs(dw{2}));   
       title('\fontsize{13} abs(dw^{(2)})');
      
elseif nDIM == 3     % Plot 3D contrast/source distribution --------------
                      % at x3 = 0 or x3 = dx/2
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 
    N3cross = floor(input.N3/2+1);
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2, Rfl{1}(:,:,N3cross));
       title('\fontsize{13} R^{(1)}=(\rho^+-\rho^-) / (\rho^++\rho^-)');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,x2,abs(dw{1}(:,:,N3cross)));
       title('\fontsize{13} abs(dw^{(1)})');
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
       IMAGESC(x1,x2+dx/2, Rfl{2}(:,:,N3cross));
       title('\fontsize{13} R^{(2)}=(\rho^+-\rho^-) / (\rho^++\rho^-)');
    subplot(1,2,2);  
       IMAGESC(x1,x2+dx/2,abs(dw{2}(:,:,N3cross)));
       title('\fontsize{13} abs(dw^{(2)})'); 
end % if 

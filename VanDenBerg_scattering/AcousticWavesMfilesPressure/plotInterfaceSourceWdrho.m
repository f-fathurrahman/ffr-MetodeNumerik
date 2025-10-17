function plotInterfaceSourceWdrho(w_drho,input)
global nDIM;  dx = input.dx;  CHI_drho = input.CHI_drho;

if nDIM == 1          % Plot 1D contrast/source distribution -------------
    x1 = input.X1; 
    set(figure,'Units','centimeters','Position',[5 5 18 12]);
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,-30:30,CHI_drho{1});  
        title('\fontsize{13} \chi^{\rho\partial\rho}_1');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,-30:30, abs(w_drho{1}));   
       title('\fontsize{13} abs(w^{\partial\rho}_1)');  
 
elseif  nDIM == 2     % Plot 2D contrast/source distribution -------------
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2,CHI_drho{1});
       title('\fontsize{13} \chi^{\rho\partial\rho}_1');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,x2, abs(w_drho{1}));   
       title('\fontsize{13} abs(w^{\partial\rho}_1)');
    set(figure,'Units','centimeters','Position',[5 5 18 12]);       
    subplot(1,2,1);  
       IMAGESC(x1,x2+dx/2,CHI_drho{2});
       title('\fontsize{13}  \chi^{\rho\partial\rho}_2');
    subplot(1,2,2);  
       IMAGESC(x1,x2+dx/2,abs(w_drho{2}));   
       title('\fontsize{13} abs(w^{\partial\rho}_2)');
      
elseif nDIM == 3     % Plot 3D contrast/source distribution --------------
                      % at x3 = 0 or x3 = dx/2
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 
    N3cross = floor(input.N3/2+1);
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2, CHI_drho{1}(:,:,N3cross));
       title('\fontsize{13} \chi^{\rho\partial\rho}_1');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,x2,abs(w_drho{1}(:,:,N3cross)));
       title('\fontsize{13} abs(w^{\partial\rho}_1)');
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
       IMAGESC(x1,x2+dx/2, CHI_drho{2}(:,:,N3cross));
       title('\fontsize{13}  \chi^{\rho\partial\rho}_2');
    subplot(1,2,2);  
       IMAGESC(x1,x2+dx/2,abs(w_drho{2}(:,:,N3cross)));
       title('\fontsize{13} abs(w^{\partial\rho}_2)'); 
       
end % if 
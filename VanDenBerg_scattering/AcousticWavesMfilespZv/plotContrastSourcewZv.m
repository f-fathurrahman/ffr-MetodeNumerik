function plotContrastSourcewZv(w_Zv,input)
global nDIM;      CHI_rho = input.CHI_rho;

if nDIM == 1           % Plot 1D contrast/source distribution -------------
    x1 = input.X1; 
    set(figure,'Units','centimeters','Position',[5 5 18 12]);
    subplot(1,2,1);  
        IMAGESC(x1,-30:30,CHI_rho);  
        title('\fontsize{13} \chi^{\rho} = 1 - \rho_{sct}/ \rho_{0}');
    subplot(1,2,2);  
        IMAGESC(x1,-30:30, abs(w_Zv{1}));   
        title('\fontsize{13} abs(w_1^{Zv})');  
 
elseif  nDIM == 2      % Plot 2D contrast/source distribution -------------
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1);  
        IMAGESC(x1,x2,abs(w_Zv{1}));
        title('\fontsize{13} abs(w_1^{Zv})');  
    subplot(1,2,2);  
        IMAGESC(x1,x2,abs(w_Zv{2}));   
        title('\fontsize{13} abs(w_2^{Zv})');  
      
elseif nDIM == 3      % Plot 3D contrast/source distribution --------------
                      % at x3 = 0 or x3 = dx/2
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 
    N3cross = floor(input.N3/2+1);
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
        IMAGESC(x1,x2,abs(w_Zv{1}(:,:,N3cross)));
        title('\fontsize{13} abs(w_1^{Zv})'); 
    subplot(1,2,2);
        IMAGESC(x1,x2,abs(w_Zv{2}(:,:,N3cross)));
        title('\fontsize{13} abs(w_2^{Zv})');        
end % if 
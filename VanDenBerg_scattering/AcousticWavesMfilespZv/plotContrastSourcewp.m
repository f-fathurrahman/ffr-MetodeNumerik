function plotContrastSourcewp(w_p,input)
global nDIM;  CHI_kap = input.CHI_kap;

if nDIM == 1           % Plot 1D contrast/source distribution -------------
    x1 = input.X1; 
    set(figure,'Units','centimeters','Position',[5 5 18 12]);
    subplot(1,2,1)
       IMAGESC(x1,-30:30, CHI_kap);  
       title('\fontsize{13} \chi^{\kappa} = 1 - \kappa_{sct}/ \kappa_{0}');
    subplot(1,2,2)
       IMAGESC(x1,-30:30, abs(w_p));   
       title('\fontsize{13} abs(w^{p})');  

elseif  nDIM == 2      % Plot 2D contrast/source distribution -------------
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1)
       IMAGESC(x1,x2, CHI_kap);  
       title('\fontsize{13} \chi^{\kappa} = 1 - \kappa_{sct}/ \kappa_{0}');
       subplot(1,2,2)
       IMAGESC(x1,x2, abs(w_p));     
       title('\fontsize{13} abs(w^p)'); 
 
elseif nDIM == 3      % Plot 3D contrast/source distribution --------------
                      % at x3 = 0 or x3 = dx/2
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 
    N3cross = floor(input.N3/2+1);
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1)   
       IMAGESC(x1,x2,CHI_kap(:,:,N3cross));
       title('\fontsize{13} \chi^{\kappa} = 1 - \kappa_{sct}/ \kappa_{0}');
    subplot(1,2,2)
       IMAGESC(x1,x2, abs(w_p(:,:,N3cross)));
       title('\fontsize{13} abs(w^p)');    
end % if 
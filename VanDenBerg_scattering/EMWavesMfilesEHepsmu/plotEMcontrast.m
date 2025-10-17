function [input] = plotEMcontrast(input)
global nDIM; 
set(figure,'Units','centimeters','Position',[5 5 18 12]); 

if  nDIM == 2;                                  % Plot 2D cross section
    x1 = input.X1(:,1);  x2 = input.X2(1,:); 
    subplot(1,2,1)   
       IMAGESC(x1,x2, input.CHI_eps);
       title('\fontsize{13}\chi^\epsilon = 1 - \epsilon_{sct}/\epsilon_0');
    subplot(1,2,2)
       IMAGESC(x1,x2, input.CHI_mu);
       title('\fontsize{13} \chi^\mu = 1 - \mu_{sct}/\mu_0');

elseif nDIM == 3;                               % Plot 3D cross section    
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 
    N3cross = floor(input.N3/2+1);         
    subplot(1,2,1)   
       IMAGESC(x1,x2, input.CHI_eps(:,:,N3cross));
       title('\fontsize{13} \chi^\epsilon = 1 - \epsilon_{sct}');
    subplot(1,2,2)
       IMAGESC(x1,x2, input.CHI_mu(:,:,N3cross));
       title('\fontsize{13} \chi^\mu = 1 - \mu_{sct}');  
       
end; % if
function [input] = plotEMcontrast(input)

set(figure,'Units','centimeters','Position',[5 5 18 12]); 
   x1 = input.X1(:,1);  x2 = input.X2(1,:); 
    subplot(1,2,1)   
       IMAGESC(x1,x2, input.CHI_eps);
       title('\fontsize{13}\chi^\epsilon = 1 - \epsilon_{sct}/\epsilon_0');
    subplot(1,2,2)
       IMAGESC(x1,x2, input.CHI_mu);
       title('\fontsize{13} \chi^\mu = 1 - \mu_{sct}/\mu_0');


function [input] = initEMcontrast(input)

input.a  = 40;      % radius circle cylinder 
       R = sqrt(input.X1.^2 + input.X2.^2);

input.CHI_eps = (1-input.eps_sct) * (R < input.a);  
input.CHI_mu  = (1-input.mu_sct)  * (R < input.a);

% Figures of contrast
x1 = input.X1(:,1);   x2 = input.X2(1,:); 
figure(1); IMAGESC(x1,x2,real(input.CHI_eps));  
title('\fontsize{13} \chi^{\epsilon} = 1 - \epsilon_{sct}/ \epsilon_{0}');
figure(2); mesh(abs(input.CHI_eps),'edgecolor', 'k'); 
           view(37.5,45); axis xy; axis('off'); axis('tight') 
           axis([1 input.N1 1 input.N2 0 1.5*max(abs(input.CHI_eps(:)))])
figure(3); IMAGESC(x1,x2,real(input.CHI_mu)); 
title('\fontsize{13}\chi^{\mu} = 1 - \mu_{sct}/ \mu_{0}');
figure(4); mesh(abs(input.CHI_mu),'edgecolor', 'k'); 
           view(37.5,45); axis xy; axis('off'); axis('tight') 
           axis([1 input.N1 1 input.N2 0 10000*max(abs(input.CHI_mu(:)))])
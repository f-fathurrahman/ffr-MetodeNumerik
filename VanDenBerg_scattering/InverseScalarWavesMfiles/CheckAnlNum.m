clear all; clc; close all; clear workspace
input = init();    s = 1 * input.dx;
   
% (1) Determine weak profile 
    P    = zeros(4);
    P(1) = input.xO(1);     P(2) = input.xO(2); 
    P(3) = input.contrast;  P(4) = input.a; 
    X1   = input.X1;        X2   = input.X2; 
    R    = sqrt((X1-P(1)).^2 + (X2-P(2)).^2);  Rp = R - P(4) + input.dx/2;
                                               Rm = R - P(4) - input.dx/2;
    CHI_weak = -P(3) * (log(1+exp(-Rp/s))-log(1+exp(-Rm/s)))/(input.dx/s);
    if sum(isnan(CHI_weak(:))) > 0 
       disp(['parameter s too small: increase value of ',num2str(s)]);
    end
    figure(13); IMAGESC(input.X1(:,1),input.X2(1,:),real(CHI_weak));  
    figure(14); mesh(abs(CHI_weak),'edgecolor', 'k'); 
                view(37.5,45); axis xy; axis('off'); axis('tight'); 
                axis([1 input.N1 1 input.N2  0 1.5*max(abs(CHI_weak(:)))])
             
% (2) Determine analytical derivatives 
    dCHI_R  =  P(3) * (1./(1+exp(Rp/s)) - 1./(1+exp(Rm/s))) / input.dx;   
    dCHI_P1 = -(X1-P(1))./R .* dCHI_R;
    dCHI_P2 = -(X2-P(2))./R .* dCHI_R;
    dCHI_P3 = CHI_weak / P(3);
    dCHI_P4 = -dCHI_R;     
  
% (3) Determine as check the numerical derivatives 
    Norm = norm(CHI_weak(:));   eps = 1e-6;      
    R = sqrt((X1-P(1)-eps).^2 + (X2-P(2)).^2); Rp = R - P(4) + input.dx/2; 
                                               Rm = R - P(4) - input.dx/2;
    dCHI_P1num = (-P(3)*(log(1+exp(-Rp/s))-log(1+exp(-Rm/s))) ...
                  /(input.dx/s) - CHI_weak) / eps;       
    disp(['Check dP1 = ',num2str( norm(dCHI_P1num(:)-dCHI_P1(:))/Norm)]); 
    R = sqrt((X1-P(1)).^2 + (X2-P(2)-eps).^2); Rp = R - P(4) + input.dx/2;
                                               Rm = R - P(4) - input.dx/2;
    dCHI_P2num = (-P(3)*(log(1+exp(-Rp/s))-log(1+exp(-Rm/s))) ...
                  /(input.dx/s) - CHI_weak) / eps; 
    disp(['Check dP2 = ',num2str( norm(dCHI_P2num(:)-dCHI_P2(:))/Norm)]); 

    R = sqrt((X1-P(1)).^2 + (X2-P(2)).^2);     Rp = R - P(4) + input.dx/2; 
                                               Rm = R - P(4) - input.dx/2;
    dCHI_P3num = (-(P(3)+eps)*(log(1+exp(-Rp/s))-log(1+exp(-Rm/s))) ...
                  /(input.dx/s) - CHI_weak) / eps;        
    disp(['Check dP3 = ',num2str( norm(dCHI_P3num(:)-dCHI_P3(:))/Norm)]); 

    R = sqrt((X1-P(1)).^2 + (X2-P(2)).^2); Rp = R - P(4)-eps + input.dx/2;
                                           Rm = R - P(4)-eps - input.dx/2;
    dCHI_P4num = (-P(3) * (log(1+exp(-Rp/s))-log(1+exp(-Rm/s))) ...
                  /(input.dx/s) - CHI_weak) / eps;
    disp(['Check dP4 = ',num2str( norm(dCHI_P4num(:)-dCHI_P4(:))/Norm)]); 
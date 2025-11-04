function [] = PSM_dct(N_mask, pz, r, pixel, k, NA, lamda, sigma, order, threshold, step, a, t_r, t_m, gamma_D, epsilon, maxloop);

% clc;
% clear;
% %%%%%%The initialization of the parameter in optical lithography system%%%%%%
% N_mask=51;   %Mask dimension
% pixel=11;   %Pixel size (nano meter)
% k=0.29;   %Process constant
% NA=1.35;   %Numerical aperture
% lamda=193;   %Wavelength (nano meter)
% sigma=0.3;   %Partial coherence factor
% order=1;   %Order of Bessel function
%threshold=18;   %Threshold to cut off the high-frequency of the DCT spectrum
                %The number of maintained low frequency coefficients is (threshold-1)*(threshold-2)/2 
                
midway=(N_mask+1)/2;   %Middle point of mask
TCC=zeros(N_mask^2,N_mask^2);   %Transmission cross coefficients
m_trinary_new=zeros(N_mask,N_mask);   %Trinary mask after the optimal post-processing

%%%%%%Calculate the transmission cross coefficients%%%%%%
disp('Calculating the transmission cross coefficients. Please wait...');
[TCC] = SOCS(N_mask, pixel, k, NA, lamda, midway, sigma, order);

%%%%%%Singular value decomposition of the partially coherent imaging system%%%%%%
[U,S,V]=svd(TCC);   %Singular value decomposition
h_1_fre=reshape(U(1:N_mask^2,1:1),N_mask,N_mask);
h_1=(fftshift(ifft2(ifftshift((h_1_fre)))));   %The impulse response of the first order approximation 
sum_eigenvalue=(sum(sum(S)));   %The summation of the eigenvalues 

%%%%%%The initialization of the parameter in the optimization%%%%%%
% step=0.2;   %Step size
% a=200;
% t_r=0.003;   %Global threshold of photoresist effect for the first eigen value
% t_m=0.33;   %Global threshold of the mask
% gamma_D=0.1;   %Weight of the discretization penalty

t_r_real=0;   %Global threshold of photoresist effect for all of the eigen value
d=zeros(N_mask,N_mask);   %Gradient of the cost function
d_D=zeros(N_mask,N_mask);   %Gradient of the discretization penalty
% epsilon=9;   %Tolerable output pattern error
% maxloop=300;   %Maximum iteration number
convergence=zeros(maxloop,1);   %Output pattern error in each iteration
count=0;   %Index of iteration number
sum6=100;   %Output pattern error corresponding to the optimized trinary mask
sum8=100;   %Output pattern error corresponding to the optimized real-valued mask

%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%
h_1_1=h_1;

for ii=1:N_mask
    for jj=1:N_mask
        h_1_vector((ii-1)*N_mask+jj)=h_1_1(ii,jj);
    end
end
for ii=1:N_mask
    for jj=1:N_mask
        g_1(ii,jj)=h_1_vector((N_mask-ii)*N_mask+(N_mask+1-jj)); %inverse vector
    end
end

%%%%%%The desired output pattern%%%%%%
% pz=zeros(N_mask,N_mask);
% for ii=16:35
%     for jj=13:23
%         pz(ii,jj)=1;
%     end
% end
% for ii=16:35
%     for jj=29:39
%         pz(ii,jj)=1;
%     end
% end
 
%%%%%%The initialization of \theta, where r=\theta%%%%%%
% r=ones(N_mask,N_mask)*pi/2;
% for ii=16:35
%     for jj=13:23
%         r(ii,jj)=pi/5;
%     end
% end
% for ii=16:35
%     for jj=29:39
%         r(ii,jj)=pi/5;
%     end
% end

%%%%%%PSM optimization in partially coherent imaging system%%%%%%
m=zeros(N_mask,N_mask);   %Mask pattern
while (sum6>epsilon) & (count<maxloop)
    count=count+1;
    %%%%%%Calculate pattern error%%%%%%
    m=cos(r);   %Gray mask
    m_trinary_p=m>t_m;
    m_trinary_n=-1*(m<(-1*t_m));
    m_trinary=m_trinary_p+m_trinary_n;   %Trinary mask
    aerial=zeros(N_mask,N_mask);   %Aerial image 
    aerial=(  abs(imfilter(double(m_trinary),h_1_1)).^2   );
    z_trinary=aerial>t_r;   %Binary output pattern
    sum6=sum(sum(abs(abs(pz)-z_trinary)));   %Output pattern error of trinary mask 
    convergence(count,1)=sum6; 
  
    %%%%%%Gradient of cost function%%%%%%
    mid1=abs(imfilter(double(m),h_1_1)).^2;
    z=1./(  1+exp(-a*mid1+a*t_r)  ); 
    mid3=(pz-z).*z.*(1-z);   
    mid4=mid3.*imfilter(double(m),h_1_1);
    mid4_4=mid3.*imfilter(double(m),conj(h_1_1));
    mid5=real(imfilter(double(mid4),conj(g_1))+imfilter(double(mid4_4),g_1));
    
    %%%%%%Gradient of discretization penaly%%%%%%  
    d_D=( (-18)*m.^3+2*m ).*((-1)*sin(r));
   
    %%%%%%%Calculate whole revision vector%%%%%%%%%%
    d=2*a*mid5.*sin(r) + gamma_D*d_D;
    r=r-step*d;   %Update
    disp(strcat('iteration=',num2str(count)));
    disp(strcat('Output pattern error = ',num2str(sum6)));
end

t_r_real=t_r*sum_eigenvalue;   %Global threshold of photoresist effect for all of the eigen value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Desired pattern%%%%%%
figure
imshow(pz,[-1,1]);
axis on;
title('Desired pattern');

%%%%%%Post-processing based on 2D-DCT%%%%%%
disp('Executing the post-processing based on 2D-DCT. Please wait...');
[m_trinary_new] = proc_dct(N_mask, pz, m, t_r, t_r_real, t_m, TCC, threshold); 

%%%%%%Save all of the data%%%%%%
save Data_PSM_dct.mat; 

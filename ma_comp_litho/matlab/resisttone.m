function [] = resisttone(N_mask, desire_pattern, distribution, pixel, k, NA, lamda, sigma, order, threshold, flag, step, a, t_r, t_m, gamma_D, gamma_WA, epsilon, maxloop);

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
% threshold=24;   %Threshold to cut off the low-frequency of the DCT spectrum
%                 %The number of maintained low frequency coefficients is
%                 %(threshold-1)*(threshold-2)/2 
% flag=0;   %The flag to control whether to include the post-processing based on 2D-DCT
%           %flag=0, the program doesn't include the post-processing
%           %flag=1, the program includeds the post-processing
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
% step=2;   %Step size
% a=200;
% t_r=0.01;   %Global threshold of photoresist effect for the first eigen value
t_r_real=0;   %Global threshold of photoresist effect for all of the eigen value
% t_m=0.33;   %Global threshold of the mask
% gamma_D=0.1;   %Weight of the discretization penalty
% gamma_WA=0.1;   %Weight of the wavelet penalty
d=zeros(N_mask,N_mask);   %Gradient of the cost function
d_D=zeros(N_mask,N_mask);   %Gradient of the discretization penalty
d_WA=zeros(N_mask,N_mask);   %Gradient of the wavelet penalty
% epsilon=8;   %Tolerable output pattern error
% maxloop=100;   %Maximum iteration number
convergence=zeros(maxloop,1);   %Output pattern error in each iteration
count=0;   %Index of iteration number
sum6=100;   %Output pattern error corresponding to the optimized trinary mask
sum8=100;   %Output pattern error corresponding to the optimized real-valued mask

%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%
h_1_1=zeros(N_mask,N_mask);
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

%%%%%%The photoresist distribution%%%%%%
% distribution=ones(N_mask,N_mask);
% for ii=16:35
%     for jj=17:18
%         distribution(ii,jj)=0;
%     end
% end
% for ii=16:35
%     for jj=34:35
%         distribution(ii,jj)=0;
%     end
% end

%%%%%%The desired pattern is a four-bars without photoresist reversion%%%%%%
% desire_pattern=zeros(N_mask,N_mask);
% for ii=16:35
%     for jj=15:16
%         desire_pattern(ii,jj)=1;
%     end
% end
% for ii=16:35
%     for jj=19:20
%         desire_pattern(ii,jj)=1;
%     end
% end
% for ii=16:35
%     for jj=32:33
%         desire_pattern(ii,jj)=1;
%     end
% end
% for ii=16:35
%     for jj=36:37
%         desire_pattern(ii,jj)=1;
%     end
% end

%%%%%%The desired pattern is two-bars with photoresist reversion%%%%%%
pz=desire_pattern-distribution+1;

%%%%%%The initialization of \theta, where r=\theta%%%%%%
r=pi/2*(pz==0) + pi/5*(pz==1);

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
    
    %%%%%%Gradient of wavelet penalty%%%%%%
    for ii=0:(N_mask/2-1)
        for jj=0:(N_mask/2-1)
            d_WA(ii*2+1,jj*2+1)= ( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+1));
            d_WA(ii*2+1,jj*2+2)= ( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+2));
            d_WA(ii*2+2,jj*2+1)= ( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+2,jj*2+1));
            d_WA(ii*2+2,jj*2+2)= ( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) * (-0.5)*sin(r(ii*2+2,jj*2+2));
        end
    end
   
    %%%%%%%Calculate whole revision vector%%%%%%%%%%
    d=2*a*mid5.*sin(r) + gamma_D*d_D + gamma_WA*d_WA;
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

if (flag==1)   %flag=0, the program doesn't include the post-processing
               %flag=1, the program includeds the post-processing
    %%%%%%Post-processing based on 2D-DCT%%%%%%
    disp('Executing the post-processing based on 2D-DCT. Please wait...');
    [m_trinary_new] = proc_dct(N_mask, pz, m, t_r, t_r_real, t_m, TCC, threshold);
end

figure
%%%%%%Photoresist distribution%%%%%%
subplot(1,3,1);
imshow(distribution,[0,1]);
title('Photoresist');
  
%%%%%%Trinary optimized mask%%%%%%
m_trinary_p=m>t_m;
m_trinary_n=-1*(m<(-1*t_m));
m_trinary=m_trinary_p+m_trinary_n;
subplot(1,3,2);
imshow(m_trinary,[-1,1]);
title('Trinary mask');
  
%%%%%%Output pattern of trinary optimized mask%%%%%%
aerial=zeros(N_mask,N_mask);
aerial_fre=zeros(N_mask,N_mask);
m_trinary_fre=(fftshift(fft2(m_trinary)));
for x=1:N_mask^2
    for y=1:N_mask^2
        index_1=mod(x-1,N_mask)+1;
        index_2=floor((x-1)/N_mask)+1;
        index_3=mod(y-1,N_mask)+1;
        index_4=floor((y-1)/N_mask)+1;
        aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)=aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)+TCC(x,y)*(m_trinary_fre(index_1,index_2))*conj(m_trinary_fre(index_3,index_4));
    end
    disp(x);
end
aerial=abs(ifft2(aerial_fre))/((N_mask)^2);
z_trinary=aerial>t_r_real;
sum6=sum(sum(abs(abs(pz)-z_trinary)));
z_trinary=z_trinary.*distribution;   %The exposed area of the negative photoresist is not etched
subplot(1,3,3);
imshow(z_trinary,[-1,1]);
xlabel(strcat('Error=',num2str(sum6)));
title('Trinary mask');

%%%%%%Convergence of optimization algorithm%%%%%%
figure   
plot([1:count],convergence(1:count,:),'k');

%%%%%%Save all of the data%%%%%%
save Data_resisttone.mat;

  
 





 
 

function [] = double_pattern(s_phi_one, s_theta_one, s_phi_two, s_theta_two, a, t_r, t_m, gamma_r_D_one, gamma_a_D_one, gamma_r_WA_one, gamma_a_WA_one, gamma_r_D_two, gamma_a_D_two, gamma_r_WA_two, gamma_a_WA_two, epsilon, maxloop);

% clc;
% clear;
%%%%%%The initialization of the parameter in the optimization%%%%%%
% s_phi_one=4;   %Step size of the \phi for the first mask 
% s_theta_one=0.01;   %Step size of the \theta for the first mask
% s_phi_two=4;   %Step size of the \phi for the second mask 
% s_theta_two=0.01;   %Step size of the \theta for the second mask
% a=80;
% t_r=0.5;   %Global threshold of the photoresist effect
% t_m=0.5;   %Global threshold of the mask amplitude
% gamma_r_D_one=0.015;   %Weight of the discretization penalty corresponding to \phi for the first mask 
% gamma_a_D_one=0.001;   %Weight of the discretization penalty corresponding to \theta for the first mask 
% gamma_r_WA_one=0.5;   %Weight of the wavelet penalty corresponding to \phi for the first mask 
% gamma_a_WA_one=0.003;   %Weight of the wavelet penalty corresponding to \theta for the first mask 
% gamma_r_D_two=0.015;   %Weight of the discretization penalty corresponding to \phi for the second mask 
% gamma_a_D_two=0.001;   %Weight of the discretization penalty corresponding to \theta for the second mask 
% gamma_r_WA_two=0.5;   %Weight of the wavelet penalty corresponding to \phi for the second mask 
% gamma_a_WA_two=0.003;   %Weight of the wavelet penalty corresponding to \theta for the second mask 
N=80;   %Mask dimension
dr_one=zeros(N,N);   %Gradient of the cost function corresponding to \phi for the first mask
da_one=zeros(N,N);   %Gradient of the cost function corresponding to \theta for the first mask
dr_D_one=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi for the first mask
da_D_one=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta for the first mask
dr_WA_one=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi for the first mask
da_WA_one=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta for the first mask
dr_two=zeros(N,N);   %Gradient of the cost function corresponding to \phi for the second mask
da_two=zeros(N,N);   %Gradient of the cost function corresponding to \theta for the second mask
dr_D_two=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi for the second mask
da_D_two=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta for the second mask
dr_WA_two=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \phi for the second mask
da_WA_two=zeros(N,N);   %Gradient of the wavelet penalty corresponding to \theta for the second mask
% epsilon=10;   %Tolerable output pattern error
% maxloop=100;   %Maximum iteration number
convergence=zeros(maxloop,1);   %Output pattern error in each iteration
count=0;   %Index of iteration number
sum6=100;   %Output pattern error corresponding to the optimized pole-level mask
cun=100;   %The lowest output pattern error that has been obtained

%%%%%%the amplitude impulse response of the coherent imaging system%%%%%%
h=fspecial('gaussian',11,14);

%%%%%%The rotation of the amplitude impulse response%%%%%%
for ii=1:11
    for j=1:11
        h1((ii-1)*11+j)=h(ii,j);
    end
end
for ii=1:11
    for j=1:11
        g(ii,j)=h1((11-ii)*11+(12-j));
    end
end

%%%%%%The desired output pattern%%%%%%
pz=zeros(N,N);
for ii=21:60
    for j=21:60
        pz(ii,j)=1;
    end
end
for ii=36:60
    for j=36:45
        pz(ii,j)=0;
    end
end

%%%%%%The initialization of \phi, where r=\phi for the first mask%%%%%%
rr_one=ones(N,N)*pi*4/5;
for ii=21:35
    for j=21:60
        rr_one(ii,j)=pi/5;
    end
end

%%%%%%The initialization of \theta, where r=\theta for the first mask%%%%%%
ra_one=ones(N,N)*pi/5;

%%%%%%The initialization of \phi, where r=\phi for the second mask%%%%%%
rr_two=ones(N,N)*pi*4/5;
for ii=36:60
    for j=21:35
        rr_two(ii,j)=pi/5;
    end
end
for ii=36:60
    for j=46:60
        rr_two(ii,j)=pi/5;
    end
end

%%%%%%The initialization of \phi, where r=\phi for the second mask%%%%%%
for ii=1:80
    for j=1:40
        ra_two(ii,j)=6*pi/5;
    end
end
for ii=1:80
    for j=41:80
        ra_two(ii,j)=pi/5;
    end
end

%%%%%%Double-patterning for generalized PSM optimization in coherent imaging system%%%%%
m_one=zeros(N,N);   %The first mask pattern
m_two=zeros(N,N);   %The second mask pattern
while (sum6>epsilon) & (count<maxloop)
    count=count+1; 
    rr_one=rr_one-s_phi_one*dr_one;   %Update
    ra_one=ra_one-s_theta_one*da_one;   %Update
    rr_two=rr_two-s_phi_two*dr_two;   %Update
    ra_two=ra_two-s_theta_two*da_two;   %Update

    m_one=0.5.*(1+cos(rr_one)).*exp(i.*ra_one);   %Calculate the first complex-valued mask pattern
    mr_one=real(m_one);   %Real part of the first complex-valued mask pattern
    mi_one=imag(m_one);   %Imaginary part of the first complex-valued mask pattern
    mmo_one=abs(m_one);   %Amplitude pattern of the first complex-valued mask pattern
    m_two=0.5.*(1+cos(rr_two)).*exp(i.*ra_two);   %Calculate the second complex-valued mask pattern
    mr_two=real(m_two);   %Real part of the second complex-valued mask pattern
    mi_two=imag(m_two);   %Imaginary part of the second complex-valued mask pattern
    mmo_two=abs(m_two);   %Amplitude pattern of the second complex-valued mask pattern  
    
    %%%%%%Quantize the first complex-valued mask to pole-level mask%%%%%%  
    viccone_one=mmo_one>t_m;
    vicctwo_one=mr_one>0;
    viccthree_one=mr_one<=0;
    viccthree_one=-1*viccthree_one;
    viccfour_one=vicctwo_one+viccthree_one;
    viccin_one=viccone_one.*viccfour_one;   %The first pole-level mask
    
    viccout_one=imfilter(viccin_one,h);
    viccbin_one=abs(viccout_one)>t_r;   %Output pattern of the first pole-level mask
    
    %%%%%%Quantize the second complex-valued mask to pole-level mask%%%%%%
    viccone_two=mmo_two>t_m;
    vicctwo_two=mr_two>0;
    viccthree_two=mr_two<=0;
    viccthree_two=-1*viccthree_two;
    viccfour_two=vicctwo_two+viccthree_two;
    viccin_two=viccone_two.*viccfour_two;   %The second pole-level mask

    viccout_two=imfilter(viccin_two,h);
    viccbin_two=abs(viccout_two)>t_r;   %Output pattern of the second pole-level mask

    sum6=sum(sum(abs(pz-((viccbin_one+viccbin_two)>=1))));   %Output pattenr error
    convergence(count,1)=sum6;
    if cun>sum6   %The lowest output pattern error that has been obtained
       cun=sum6;
    end
    disp(cun);

    mid1_one=imfilter(m_one,h);   %Convolution between the first complex-valued mask and low-pass filter
    mid1mo_one=abs(mid1_one);   %Convolution between the first complex-valued mask amplitude and low-pass filter
    mid1r_one=imfilter(mr_one,h);   %Convolution between real part of the first complex-valued mask amplitude and low-pass filter
    mid1i_one=imfilter(mi_one,h);   %Convolution between imaginary part of the first complex-valued mask amplitude and low-pass filter
    mid1_two=imfilter(m_two,h);   %Convolution between the second complex-valued mask and low-pass filter
    mid1mo_two=abs(mid1_two);   %Convolution between the second complex-valued mask amplitude and low-pass filter
    mid1r_two=imfilter(mr_two,h);   %Convolution between real part of the second complex-valued mask amplitude and low-pass filter
    mid1i_two=imfilter(mi_two,h);   %Convolution between imaginary part of the second complex-valued mask amplitude and low-pass filter
   
    z_one=1./ (  1+exp(-1*a*(mid1mo_one)+a*t_r)  ); 
    z_two=1./ (  1+exp(-1*a*(mid1mo_two)+a*t_r)  ); 

    F=0.5*(tanh(z_one+z_two-1)+1);   %Overall photoresist effect of the two exposures

    %%%%%%Calculations for the gradient of the first mask%%%%%%
    mid3_one=(pz-F).*(sech(z_one+z_two-1)).^2.*z_one.*(1-z_one).*mid1r_one.*(1./mid1mo_one);   
    mid5_one=imfilter(mid3_one,g);  
    mid7_one=(pz-F).*(sech(z_one+z_two-1)).^2.*z_one.*(1-z_one).*mid1i_one.*(1./mid1mo_one);   
    mid9_one=imfilter(mid7_one,g);
    %%%%%%Calculations for the gradient of the second mask%%%%%%
    mid3_two=(pz-F).*(sech(z_one+z_two-1)).^2.*z_two.*(1-z_two).*mid1r_two.*(1./mid1mo_two);   
    mid5_two=imfilter(mid3_two,g);   
    mid7_two=(pz-F).*(sech(z_one+z_two-1)).^2.*z_two.*(1-z_two).*mid1i_two.*(1./mid1mo_two);   
    mid9_two=imfilter(mid7_two,g);
    
    %%%%%%Gradient of discretization penaly of the first mask%%%%%%
    dr_D_one=(-0.5)*sin(rr_one).*(1+cos(rr_one));
    da_D_one=8*( sin(4*ra_one-pi*3/2)+1 ).*cos(4*ra_one-pi*3/2);
 
    %%%%%%Gradient of discretization penaly of the second mask%%%%%%
    dr_D_two=(-0.5)*sin(rr_two).*(1+cos(rr_two));
    da_D_two=8*( sin(4*ra_two-pi*3/2)+1 ).*cos(4*ra_two-pi*3/2);
    
    %%%%%%Gradient of wavelet penalty of the first mask%%%%%%
    for ii=0:N/2-1
        for jj=0:N/2-1
            dr_WA_one(ii*2+1,jj*2+1)=  -1*sin(rr_one(ii*2+1,jj*2+1))*real( exp((-i)*ra_one(ii*2+1,jj*2+1)) *( 3*m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            dr_WA_one(ii*2+1,jj*2+2)=  -1*sin(rr_one(ii*2+1,jj*2+2))*real( exp((-i)*ra_one(ii*2+1,jj*2+2)) *( 3*m_one(ii*2+1,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            dr_WA_one(ii*2+2,jj*2+1)=  -1*sin(rr_one(ii*2+2,jj*2+1))*real( exp((-i)*ra_one(ii*2+2,jj*2+1)) *( 3*m_one(ii*2+2,jj*2+1) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+2) ) );
            dr_WA_one(ii*2+2,jj*2+2)=  -1*sin(rr_one(ii*2+2,jj*2+2))*real( exp((-i)*ra_one(ii*2+2,jj*2+2)) *( 3*m_one(ii*2+2,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) ) );
        end
    end
    for ii=0:N/2-1
        for jj=0:N/2-1
            da_WA_one(ii*2+1,jj*2+1)=  (1+cos(rr_one(ii*2+1,jj*2+1)))*real( (-i)*exp((-i)*ra_one(ii*2+1,jj*2+1)) *( 3*m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            da_WA_one(ii*2+1,jj*2+2)=  (1+cos(rr_one(ii*2+1,jj*2+2)))*real( (-i)*exp((-i)*ra_one(ii*2+1,jj*2+2)) *( 3*m_one(ii*2+1,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+2,jj*2+1) - m_one(ii*2+2,jj*2+2) ) );
            da_WA_one(ii*2+2,jj*2+1)=  (1+cos(rr_one(ii*2+2,jj*2+1)))*real( (-i)*exp((-i)*ra_one(ii*2+2,jj*2+1)) *( 3*m_one(ii*2+2,jj*2+1) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+2) ) );
            da_WA_one(ii*2+2,jj*2+2)=  (1+cos(rr_one(ii*2+2,jj*2+2)))*real( (-i)*exp((-i)*ra_one(ii*2+2,jj*2+2)) *( 3*m_one(ii*2+2,jj*2+2) - m_one(ii*2+1,jj*2+1) - m_one(ii*2+1,jj*2+2) - m_one(ii*2+2,jj*2+1) ) );
        end
    end
    
    %%%%%%Gradient of wavelet penalty of the second mask%%%%%%
    for ii=0:N/2-1
        for jj=0:N/2-1
            dr_WA_two(ii*2+1,jj*2+1)=  -1*sin(rr_two(ii*2+1,jj*2+1))*real( exp((-i)*ra_two(ii*2+1,jj*2+1)) *( 3*m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            dr_WA_two(ii*2+1,jj*2+2)=  -1*sin(rr_two(ii*2+1,jj*2+2))*real( exp((-i)*ra_two(ii*2+1,jj*2+2)) *( 3*m_two(ii*2+1,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            dr_WA_two(ii*2+2,jj*2+1)=  -1*sin(rr_two(ii*2+2,jj*2+1))*real( exp((-i)*ra_two(ii*2+2,jj*2+1)) *( 3*m_two(ii*2+2,jj*2+1) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+2) ) );
            dr_WA_two(ii*2+2,jj*2+2)=  -1*sin(rr_two(ii*2+2,jj*2+2))*real( exp((-i)*ra_two(ii*2+2,jj*2+2)) *( 3*m_two(ii*2+2,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) ) );
        end
    end
    for ii=0:N/2-1
        for jj=0:N/2-1
            da_WA_two(ii*2+1,jj*2+1)=  (1+cos(rr_two(ii*2+1,jj*2+1)))*real( (-i)*exp((-i)*ra_two(ii*2+1,jj*2+1)) *( 3*m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            da_WA_two(ii*2+1,jj*2+2)=  (1+cos(rr_two(ii*2+1,jj*2+2)))*real( (-i)*exp((-i)*ra_two(ii*2+1,jj*2+2)) *( 3*m_two(ii*2+1,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+2,jj*2+1) - m_two(ii*2+2,jj*2+2) ) );
            da_WA_two(ii*2+2,jj*2+1)=  (1+cos(rr_two(ii*2+2,jj*2+1)))*real( (-i)*exp((-i)*ra_two(ii*2+2,jj*2+1)) *( 3*m_two(ii*2+2,jj*2+1) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+2) ) );
            da_WA_two(ii*2+2,jj*2+2)=  (1+cos(rr_two(ii*2+2,jj*2+2)))*real( (-i)*exp((-i)*ra_two(ii*2+2,jj*2+2)) *( 3*m_two(ii*2+2,jj*2+2) - m_two(ii*2+1,jj*2+1) - m_two(ii*2+1,jj*2+2) - m_two(ii*2+2,jj*2+1) ) );
        end
    end

    %%%%%%Gradient of the overall cost function for the first mask%%%%%%
    dr_one=0.5*( a*sin(rr_one).*cos(ra_one).*mid5_one + a*sin(rr_one).*sin(ra_one).*mid9_one ) + gamma_r_D_one*dr_D_one + gamma_r_WA_one*dr_WA_one;
    da_one=a*0.5*(1+cos(rr_one)).*sin(ra_one).*mid5_one - a*0.5*(1+cos(rr_one)).*cos(ra_one).*mid9_one + gamma_a_D_one*da_D_one + gamma_a_WA_one*da_WA_one;      

    %%%%%%Gradient of the overall cost function for the second mask%%%%%%
    dr_two=0.5*( a*sin(rr_two).*cos(ra_two).*mid5_two + a*sin(rr_two).*sin(ra_two).*mid9_two ) + gamma_r_D_two*dr_D_two + gamma_r_WA_two*dr_WA_two;
    da_two=a*0.5*(1+cos(rr_two)).*sin(ra_two).*mid5_two - a*0.5*(1+cos(rr_two)).*cos(ra_two).*mid9_two + gamma_a_D_two*da_D_two + gamma_a_WA_two*da_WA_two; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Desired pattern%%%%%%
figure
imshow(pz,[-1,1]);
colormap('default');
axis on;
title('Desired pattern');

%%%%%%The first mask%%%%%%
figure
subplot(1,3,1);
imshow(viccin_one,[-1,1]);
colormap('default');
axis on;
title('The first mask');

%%%%%%The second mask%%%%%%
subplot(1,3,2);
imshow(viccin_two,[-1,1]);
colormap('default');
axis on;
title('The second mask');

%%%%%%The output pattern%%%%%%
subplot(1,3,3);
imshow(double((viccbin_one+viccbin_two)>=1),[-1,1]);
colormap('default');
axis on;
title('The output pattern');
xlabel(strcat('Error=',num2str(sum6)));

%%%%%%Convergence of optimization algorithm%%%%%%
figure   
plot([1:count],convergence(1:count,:),'k');

%%%%%%Save all of the data%%%%%%
save Data_double_pattern.mat;
        

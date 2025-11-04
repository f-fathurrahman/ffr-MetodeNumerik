function [] = GPSM_tv(N, pz, ra, phase_n, s_phi, s_theta, a, t_r, t_m, gamma_r_D, gamma_a_D, gamma_r_TV, gamma_a_TV, epsilon, maxloop);

% clc;
% clear;
% %%%%%%The initialization of the parameter in the optimization%%%%%%
% N=80;   %Mask dimension
% s_phi=2;   %Step size of the \phi
% s_theta=0.01;   %Step size of the \theta
% a=80;
% t_r=0.5;   %Global threshold of the photoresist effect
% t_m=0.5;   %Global threshold of the mask amplitude
% gamma_r_D=0.01;   %Weight of the discretization penalty corresponding to \phi
% gamma_a_D=0.001;   %Weight of the discretization penalty corresponding to \theta
% gamma_r_TV=0.1;   %Weight of the total variation penalty corresponding to \phi
% gamma_a_TV=0.001;   %Weight of the total variation penalty corresponding to \theta
dr=zeros(N,N);   %Gradient of the cost function corresponding to \phi
da=zeros(N,N);   %Gradient of the cost function corresponding to \theta
dr_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \phi
da_D=zeros(N,N);   %Gradient of the discretization penalty corresponding to \theta
dr_TV=zeros(N,N);   %Gradient of the total variation penalty corresponding to \phi
da_TV=zeros(N,N);   %Gradient of the total variation penalty corresponding to \theta
% epsilon=44;   %Tolerable output pattern error
% maxloop=27;   %Maximum iteration number
convergence=zeros(maxloop,1);   %Output pattern error in each iteration
count=0;   %Index of iteration number
sum6=100;   %Output pattern error corresponding to the optimized pole-level mask
sum8=100;   %Output pattern error corresponding to the optimized complex-valued mask
depth=1;
wavename='haar';

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
% pz=zeros(N,N);
% for ii=21:60
%     for j=21:60
%         pz(ii,j)=1;
%     end
% end
% for ii=36:60
%     for j=36:45
%         pz(ii,j)=0;
%     end
% end
figure
imshow(pz,[-1,1]);
colormap('default');
axis on;

%%%%%%The initialization of \phi, where r=\phi%%%%%%
rr=pi*4/5*(pz==0) + pi/5*(pz==1);

% %%%%%%The initialization of \theta, where r=\theta%%%%%%
% for ii=1:35
%     for j=36:45
%         ra(ii,j)=7*pi/4;
%     end
% end
% for ii=21:35
%     for j=36:45
%         ra(ii,j)=pi/4;
%     end
% end
% %%%%%%top left%%%%%%
% for ii=1:35
%     for j=1:35
%         ra(ii,j)=7*pi/4;
%     end
% end
% for ii=21:35
%     for j=21:35
%         ra(ii,j)=pi/4;
%     end
% end
% %%%%%%top right%%%%%%
% for ii=1:35
%     for j=46:80
%         ra(ii,j)=7*pi/4;
%     end
% end
% for ii=21:35
%     for j=46:60
%         ra(ii,j)=pi/4;
%     end
% end
% %%%%%%left%%%%%%
% for ii=36:80
%     for j=1:40
%         ra(ii,j)=5*pi/4;
%     end
% end
% for ii=36:60
%     for j=21:35
%         ra(ii,j)=3*pi/4;
%     end
% end
% %%%%%%right%%%%%%
% for ii=36:80
%     for j=41:80
%         ra(ii,j)=pi/4;
%     end
% end
% for ii=36:60
%     for j=46:60
%         ra(ii,j)=7*pi/4;
%     end
% end

%%%%%%Generalized PSM optimization in coherent imaging system%%%%%
m=zeros(N,N);   %Mask pattern
cun=1000;
while (sum6>epsilon) & (count<maxloop)
    count=count+1; 
    rr=rr-s_phi*dr;   %Update
    ra=ra-s_theta*da;   %Update
    m=0.5.*(1+cos(rr)).*exp(i.*ra);   %Calculate complex-valued mask pattern
    mr=real(m);   %Real part of complex-valued mask pattern
    mi=imag(m);   %Imaginary part of complex-valued mask pattern
    mmo=abs(m);   %Amplitude pattern of complex-valued mask pattern
    
    %%%%%%Quantize the complex-valued mask to pole-level mask%%%%%
    if (phase_n==4)   %Four-phase PSM
        viccone=mmo>t_m;   %Transparent area on the mask
        vicctwo=(mr>=0)&(mi>=0);   %Area with phase of pi/4
        vicctwo=exp(i*pi/4)*vicctwo;
        viccthree=(mr<=0)&(mi>=0);   %Area with phase of pi*3/4
        viccthree=exp(i*pi*3/4)*viccthree;
        viccfour=(mr<=0)&(mi<=0);   %Area with phase of pi*5/4
        viccfour=exp(i*pi*5/4)*viccfour;
        viccfive=(mr>=0)&(mi<=0);   %Area with phase of pi*7/4
        viccfive=exp(i*pi*7/4)*viccfive;
        viccsix=vicctwo+viccthree+viccfour+viccfive;   %Phase pattern of mask pattern
        viccin=viccone.*viccsix;   %Pole-level mask pattern
    elseif (phase_n==2)   %Two-phase PSM
        viccone=mmo>t_m;   %Transparent area on the mask
        vicctwo=mr>0;   %Area with phase of 0
        viccthree=mr<=0;   %Area with phase of pi
        viccthree=-1*viccthree;
        viccfour=vicctwo+viccthree;   %Phase pattern of mask pattern
        viccin=viccone.*viccfour;   %Pole-level mask pattern
    end
    
    viccout=imfilter(viccin,h);
    viccbin=abs(viccout)>t_r;   %Output pattern of pole-level mask

    sum6=sum(sum(abs(pz-viccbin)));    
    convergence(count,1)=sum6;
    
    if cun>sum6
        cun=sum6;
    end
    disp(cun);
   
    mid1=imfilter(m,h);   %Convolution between complex-valued mask and low-pass filter
    mid1mo=abs(mid1);   %Convolution between complex-valued mask amplitude and low-pass filter
    mid1r=imfilter(mr,h);   %Convolution between real part of complex-valued mask amplitude and low-pass filter 
    mid1i=imfilter(mi,h);   %Convolution between imaginary part of complex-valued mask amplitude and low-pass filter
    z=1./ (  1+exp(-1*a*(mid1mo)+a*t_r)  ); 
    mid3=( pz-z ).*z.*(1-z).*mid1r.*(1./mid1mo);   
    mid5=imfilter(mid3,g);   
    mid7=( pz-z ).*z.*(1-z).*mid1i.*(1./mid1mo);   
    mid9=imfilter(mid7,g);
    
    %%%%%%Gradient of the discretization penalty corresponding to \phi%%%%%%  
    dr_D=(-0.5)*sin(rr).*(1+cos(rr));
    
    %%%%%%Gradient of the discretization penalty corresponding to \theta%%%%%%  
    if (phase_n==4)   %Four-phase PSM
        da_D=8*( sin(4*ra-pi*3/2) + 1 ).*cos(4*ra-pi*3/2);
    elseif (phase_n==2)   %Two-phase PSM
        da_D=4.*( sin(2.*ra-pi/2)+1 ).*cos(2.*ra-pi/2);
    end
    
    %%%%%%Gradient of total variation penaly%%%%%%
    f=abs(m-pz);
    f_right=zeros(N,N);   %Right shift of f
    f_right(1:N,2:N)=f(1:N,1:N-1);
    f_up=zeros(N,N);   %Up shift of f
    f_up(1:N-1,1:N)=f(2:N,1:N);
    f1=sign( f-f_right );
    f1(:,1)=0;
    f2=sign( f-f_up );
    f2(N,:)=0;
   
    f1_left=zeros(N,N);   %Left shift of f1
    f1_left(1:N,1:N-1)=f1(1:N,2:N);
    f2_down=zeros(N,N);   %Down shift of f2
    f2_down(2:N,N:N)=f2(1:N-1,N:N);
    f11=f1-f1_left;
    f11(:,N)=0;
    f22=f2-f2_down;
    f22(1,:)=0;    
    
    f3=(f11+f22);
    %%%%%%Gradient of total variation penaly corresponding to \phi%%%%%%
    dr_TV=f3*(-1).*sin(rr).*real( (m-z).*exp((-1)*i*ra) ) /2 ./abs(m-z);
    %%%%%%Gradient of total variation penaly corresponding to \theta%%%%%%
    da_TV=f3.*(1+cos(rr)).*real( (m-z).*exp((-1)*i*ra)*(-1)*i ) /2 ./abs(m-z);

    %%%%%%Gradient of overall cost function%%%%%% 
    dr=a*sin(rr).*cos(ra).*mid5 + a*sin(rr).*sin(ra).*mid9 + gamma_r_D*dr_D + gamma_r_TV*dr_TV;
    da=2*a*0.5*(1+cos(rr)).*sin(ra).*mid5 - 2*a*0.5*(1+cos(rr)).*cos(ra).*mid9 + gamma_a_D*da_D + gamma_a_TV*da_TV;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Output pattern of desired pattern%%%%%%
output_pz=imfilter(pz,h)>t_r;
figure   %Output pattern of desired pattern
imshow(output_pz,[-1,1]);
colormap('default');
axis on;
sum10=sum(sum(abs(pz-output_pz)));   %Output pattern error of desired pattern

%%%%%%Magnitude of optimized complex-valued mask%%%%%%
figure
imshow((cos(rr)+1)/2,[-1,1]);
colormap('default');
axis on;
colorbar('vertical');

%%%%%%Phase of optimized complex-valued mask%%%%%%
majorangle=mod(ra,2*pi);   %The phase restricted in [0,2pi]
figure
imshow(majorangle,[0,2*pi]);
colormap('default');
axis on;
colorbar('vertical');

%%%%%%Output of optimized complex-valued mask%%%%%%
output_m=abs(imfilter(m,h))>t_r;
figure   
imshow(output_m,[-1,1]);
colormap('default');
axis on;
sum8=sum(sum(abs(pz-output_m)));   %Output pattern error of optimized complex-valued mask

%%%%%%Magnitude of optimized pole-level mask%%%%%%
m=0.5.*(1+cos(rr)).*exp(i.*ra);   %Calculate complex-valued mask pattern
mr=real(m);   %Real part of complex-valued mask pattern
mi=imag(m);   %Imaginary part of complex-valued mask pattern
mmo=abs(m);   %Amplitude pattern of complex-valued mask pattern
viccone=mmo>t_m;   %Transparent area on the mask
figure
imshow(abs(viccone),[-1,1]);
colormap('default');
axis on;

%%%%%%Phase of optimized pole-level mask%%%%%%
if (phase_n==4)   %Four-phase PSM
    vicconeangle=mmo>t_m;   %Transparent area on the mask
    vicctwoangle=(mr>=0)&(mi>=0);   %Area with phase of pi/4
    vicctwoangle=(pi/4)*vicctwoangle;
    viccthreeangle=(mr<=0)&(mi>=0);   %Area with phase of pi*3/4
    viccthreeangle=(pi*3/4)*viccthreeangle;
    viccfourangle=(mr<=0)&(mi<=0);   %Area with phase of pi*5/4
    viccfourangle=(pi*5/4)*viccfourangle;
    viccfiveangle=(mr>=0)&(mi<=0);   %Area with phase of pi*7/4
    viccfiveangle=(pi*7/4)*viccfiveangle;
    viccsixangle=vicctwoangle+viccthreeangle+viccfourangle+viccfiveangle;   %Phase pattern of mask pattern
    viccinangle=vicconeangle.*viccsixangle;
elseif (phase_n==2)   %Two-phase PSM
    vicconeangle=mmo>t_m;   %Transparent area on the mask
    vicctwoangle=mr>0;   %Area with phase of 0
    vicctwoangle=0*vicctwoangle;
    viccthreeangle=mr<=0;   %Area with phase of pi
    viccthreeangle=pi*viccthreeangle;
    viccfourangle=vicctwoangle+viccthreeangle;   %Phase pattern of mask pattern
    viccinangle=vicconeangle.*viccfourangle;
end
figure
imshow(viccinangle,[0,2*pi]);
colormap('default');
axis on;
colorbar('vertical');

%%%%%%Output of optimized pole-level mask%%%%%%
if (phase_n==4)   %Four-phase PSM
    viccone=mmo>t_m;   %Transparent area on the mask
    vicctwo=(mr>=0)&(mi>=0);   %Area with phase of pi/4
    vicctwo=exp(i*pi/4)*vicctwo;
    viccthree=(mr<=0)&(mi>=0);   %Area with phase of pi*3/4
    viccthree=exp(i*pi*3/4)*viccthree;
    viccfour=(mr<=0)&(mi<=0);   %Area with phase of pi*5/4
    viccfour=exp(i*pi*5/4)*viccfour;
    viccfive=(mr>=0)&(mi<=0);   %Area with phase of pi*7/4
    viccfive=exp(i*pi*7/4)*viccfive;
    viccsix=vicctwo+viccthree+viccfour+viccfive;   %Phase pattern of mask pattern
    viccin=viccone.*viccsix;   %Pole-level mask pattern
elseif (phase_n==2)   %Four-phase PSM
    viccone=mmo>t_m;   %Transparent area on the mask
    vicctwo=mr>0;   %Area with phase of 0
    viccthree=mr<=0;   %Area with phase of pi
    viccthree=-1*viccthree;
    viccfour=vicctwo+viccthree;   %Phase pattern of mask pattern
    viccin=viccone.*viccfour;   %Pole-level mask pattern
end

output_m_b=abs(imfilter(viccin,h))>t_r;
figure   %Output pattern of optimized pole-level mask 
imshow(output_m_b,[-1,1]);
colormap('default');
axis on;
sum6=sum(sum(abs(pz-output_m_b)));   %Output pattern error of optimized pole-level mask

%%%%%%Convergence of optimization algorithm%%%%%%
figure   
plot([1:count],convergence(1:count,:),'k');


%%%%%%Wavelet transform of amplitude of optimized complex-valued mask%%%%%%
[coefs,sizes] = wavedec2(abs(viccin)*2-1,depth,wavename);
coef_a = appcoef2(coefs,sizes,wavename,depth);
[coef_h, coef_v, coef_d] = detcoef2('all',coefs,sizes,depth);
waveletimage=zeros(80,80);
waveletimage(1:N/2,1:N/2)=coef_a;
waveletimage(1:N/2,N/2+1:N)=coef_h;
waveletimage(N/2+1:N,1:N/2)=coef_v;
waveletimage(N/2+1:N,N/2+1:N)=coef_d;
figure
subplot(1,2,1)
imshow(abs(viccin),[-1,1]);
colormap('default');
axis on;
colorbar('horizantal');
subplot(1,2,2)
imshow(waveletimage,[min(min(waveletimage)),max(max(waveletimage))]);
colormap('default');
axis on;
colorbar('horizantal');

%%%%%%Save all of the data%%%%%%
save Data_GPSM_tv.mat;



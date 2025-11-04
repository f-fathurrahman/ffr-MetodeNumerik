function [] = OPC_fse(N, pz, N_filter, pixel, k, NA, lamda, order, sigma_large_inner, sigma_large_outer, step, a, t_r, t_m, gamma_D, gamma_WA, epsilon, maxloop);

% clc;
% clear;
% %%%%%%The initialization of the parameter in the optimization%%%%%%
% N=184;   %Mask dimension
% N_filter=21;   %Amplitude impulse response dimension
% pixel=5.625;   %Pixel size (nano meter)
% k=0.29;
% NA=1.25;
% lamda=193;   %Wavelength (nano meter)
% order=1;   %Order of Bessel function
% sigma_large_inner=0.8;   %Inner Partial coherence factor
% sigma_large_outer=0.975;   %Outer Partial coherence factor

D=pixel*184;
D_C_1=lamda/2/sigma_large_outer/NA;   %coherence length
D_C_2=lamda/2/sigma_large_inner/NA;   %coherence length
omega_0=pi/D;
N_coherence=floor(2*D/(2*D_C_1)+1)+2;   %Size of annular illumination
midway_coherence=(N_coherence+1)/2;   %Middle point of annular illumination

% step=2;   %Step size
% a=25;
% t_r=0.19;
% t_m=0.5;   %Global threshold of the mask
% gamma_D=0.025;   %Weight of the discretization penalty
% gamma_WA=0.025;   %Weight of the wavelet penalty
d=zeros(N,N);   %Gradient of the cost function
d_D=zeros(N,N);   %Gradient of the discretization penalty
d_WA=zeros(N,N);   %Gradient of the wavelet penalty
% epsilon=0;   %Tolerable output pattern error
% maxloop=43;   %Maximum iteration number
convergence=zeros(maxloop,1);   %Output pattern error in each iteration
count=0;   %Index of iteration number
sum6=10000;   %Output pattern error corresponding to the optimized binary mask
sum8=10000;   %Output pattern error corresponding to the optimized real-valued mask


%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%
h=zeros(N_filter,N_filter);
radius=0;
midway=(N_filter+1)/2;%middle point of low pass filter
for row=1:N_filter
    for column=1:N_filter
        radius=pixel*sqrt( (row-midway)^2 + (column-midway)^2 );
        if (radius<=(midway)*pixel)
            argument=2*pi*radius*NA/lamda;
            if (radius==0)
                h(row,column)=h(row-1,column);
            else
                h(row,column)=besselj(order,argument)/argument;
            end
        end
    end
end
h=h/sum(sum(h));   %Normalization
figure
surf(h);
figure
surf(abs(fftshift(fft2(h))));

for ii=1:N_filter
    for jj=1:N_filter
        h_vector((ii-1)*N_filter+jj)=h(ii,jj);
    end
end
for ii=1:N_filter
    for jj=1:N_filter
        g(ii,jj)=h_vector((N_filter-ii)*N_filter+(N_filter+1-jj)); %inverse vector
    end
end
figure
surf(g);
 
%%%%%%The desired output pattern%%%%%%
% pz=zeros(N,N);
% for ii=9:16
%     for jj=11:174
%         pz(ii,jj)=1;
%     end
% end
% for ii=25:32
%     for jj=31:70
%         pz(ii,jj)=1;
%     end
% end
% for ii=33:40
%     for jj=31:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=41:48
%     for jj=11:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=73:80
%     for jj=31:70
%         pz(ii,jj)=1;
%     end
% end
% for ii=65:72
%     for jj=31:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=57:64
%     for jj=11:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=25:32
%     for jj=135:174
%         pz(ii,jj)=1;
%     end
% end
% for ii=33:40
%     for jj=135:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=41:48
%     for jj=115:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=73:80
%     for jj=135:174
%         pz(ii,jj)=1;
%     end
% end
% for ii=65:72
%     for jj=135:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=57:64
%     for jj=115:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=89:96
%     for jj=11:174
%         pz(ii,jj)=1;
%     end
% end
% for ii=105:112
%     for jj=31:70
%         pz(ii,jj)=1;
%     end
% end
% for ii=113:120
%     for jj=31:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=121:128
%     for jj=11:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=153:160
%     for jj=31:70
%         pz(ii,jj)=1;
%     end
% end
% for ii=145:152
%     for jj=31:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=137:144
%     for jj=11:50
%         pz(ii,jj)=1;
%     end
% end
% for ii=105:112
%     for jj=135:174
%         pz(ii,jj)=1;
%     end
% end
% for ii=113:120
%     for jj=135:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=121:128
%     for jj=115:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=153:160
%     for jj=135:174
%         pz(ii,jj)=1;
%     end
% end
% for ii=145:152
%     for jj=135:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=137:144
%     for jj=115:154
%         pz(ii,jj)=1;
%     end
% end
% for ii=169:176
%     for jj=11:174
%         pz(ii,jj)=1;
%     end
% end
figure
imshow(pz);
 
%%%%%%The initialization of \theta, where r=\theta%%%%%%
r=pi*4/5*(pz==0) + pi/5*(pz==1);
 
%%%%%%%Annular illumination pattern%%%%%%%
radius_1=D/(2*D_C_1);   %Inner radius of annular illumination
radius_2=D/(2*D_C_2);   %Outer radius of annular illumination 
yita=zeros(N_coherence,N_coherence);   %Annular illumination pattern
for row=1:N_coherence
    for column=1:N_coherence
        radius=pixel*sqrt( (row-midway_coherence)^2 + (column-midway_coherence)^2 );
        if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
            yita(row,column)=1;
        end
    end
end
yita=yita/sum(sum(yita));
figure
imshow(yita>0);

%%%%%%OPC optimization in partially coherent imaging system%%%%%%
m=zeros(N,N);   %Mask pattern
while (sum6>epsilon) & (count<maxloop)   
   count=count+1; 
   r=r-step*d;   %Update
   %%%%%%Calculate pattern error%%%%%%
   m=(1+cos(r))/2;   %Grey mask
   m_binary=m>t_m;   %Binary mask
   aerial=zeros(N,N);   %Aerial image
   for p=1:N_coherence
       for q=1:N_coherence
           radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
           if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
               exponential=zeros(N_filter,N_filter); % has the same dimension as the filter h
               for row=1:N_filter
                   for column=1:N_filter
                       argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                       exponential(row,column)=exp(i*omega_0*argument);
                   end
               end
               aerial=aerial+yita(p,q)* abs(  imfilter(double(m_binary),h.*exponential)  ).^2;
           end
       end
   end
   z_binary=aerial>t_r;   %Binary output pattern
   sum6=sum(sum(abs(pz-double(z_binary))));
   convergence(count,1)=sum6; 
   
   %%%%%%Gradient of cost function%%%%%%
   mid1=zeros(N,N);
   for p=1:N_coherence
       for q=1:N_coherence
           radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
           if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
               exponential=zeros(N_filter,N_filter);
               for row=1:N_filter
                   for column=1:N_filter
                       argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                       exponential(row,column)=exp(i*omega_0*argument);
                   end
               end
               mid1=mid1+yita(p,q)* abs(imfilter(double(m),h.*exponential)).^2;
           end
       end
   end
   z=1./ (  1+exp(-a*mid1+a*t_r)  ); 
   mid3=( pz-z ).*z.*(1-z);   
   mid5=zeros(N,N);
   for p=1:N_coherence
       for q=1:N_coherence
           radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
           if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
               exponential=zeros(N_filter,N_filter);
               for row=1:N_filter
                   for column=1:N_filter
                       argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                       exponential(row,column)=exp(i*omega_0*argument);
                   end
               end
               mid5=mid5+yita(p,q)* real(  imfilter( mid3.*imfilter(double(m),h.*exponential) , conj(g.*exponential) )  +   imfilter( mid3.*imfilter(double(m),conj(h.*exponential)) , g.*exponential )    );
           end
       end
   end
   
   %%%%%%Gradient of discretization penaly%%%%%%  
   d_D=( (-8)*m+4 )*(-0.5).*sin(r);
   
   %%%%%%Gradient of wavelet penaly%%%%%%
   for ii=0:(N/2-1)
       for jj=0:(N/2-1)
           d_WA(ii*2+1,jj*2+1)= ( 3*m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+1));
           d_WA(ii*2+1,jj*2+2)= ( 3*m(ii*2+1,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+2,jj*2+1) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+1,jj*2+2));
           d_WA(ii*2+2,jj*2+1)= ( 3*m(ii*2+2,jj*2+1) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+2) ) * (-0.5)*sin(r(ii*2+2,jj*2+1));
           d_WA(ii*2+2,jj*2+2)= ( 3*m(ii*2+2,jj*2+2) - m(ii*2+1,jj*2+1) - m(ii*2+1,jj*2+2) - m(ii*2+2,jj*2+1) ) * (-0.5)*sin(r(ii*2+2,jj*2+2));
       end
   end
 
   %%%%%%Gradient of overall cost function%%%%%%
   d=a*mid5.*(0.5*sin(r))+gamma_D*d_D +gamma_WA*d_WA;

   disp(count);
   disp(sum6);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Desired pattern%%%%%%
figure
subplot(2,3,1);
imshow(pz);
title('Desired pattern');

%%%%%%Output pattern of desired pattern%%%%%%
mid=zeros(N,N);
for p=1:N_coherence
    for q=1:N_coherence
        radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
        if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
            exponential=zeros(N_filter,N_filter);
            for row=1:N_filter
                for column=1:N_filter
                    argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                    exponential(row,column)=exp(i*omega_0*argument);
                end
            end
            mid=mid+yita(p,q)* abs(  imfilter(double(pz),h.*exponential)  ).^2;
        end
    end
end
z4=mid >t_r;
sum10=sum(sum(abs(pz-double(z4))));   %Output pattern error of desired pattern
subplot(2,3,4);
imshow(z4);
xlabel(strcat('error=',num2str(sum10)));

%%%%%%Real-valued optimized mask%%%%%%
subplot(2,3,2);
imshow(m);
title('Real-valued optimized mask');

%%%%%%Output pattern of real-valued optimized mask%%%%%%
grayout=zeros(N,N);
for p=1:N_coherence
    for q=1:N_coherence
        radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
        if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
            exponential=zeros(N_filter,N_filter);
            for row=1:N_filter
                for column=1:N_filter
                    argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                    exponential(row,column)=exp(i*omega_0*argument);
                end
            end
            grayout=grayout+yita(p,q)* abs(  imfilter(double(m),h.*exponential)  ).^2;
        end
    end
end
binout=grayout>t_r;
sum8=sum(sum(abs(pz-double(binout))));   %Output pattern error of real-valued optimized mask
subplot(2,3,5);
imshow(binout);
xlabel(strcat('error=',num2str(sum8)));

%%%%%%Binary optimized mask%%%%%%
m_binary=m>t_m;
subplot(2,3,3);
imshow(m_binary);
title('Binary optimized mask');

%%%%%%Output pattern of binary optimized mask%%%%%%
aerial_1=zeros(N,N);
for p=1:N_coherence
    for q=1:N_coherence 
        radius=pixel*sqrt( (p-midway_coherence)^2 + (q-midway_coherence)^2 );    
        if (radius<=radius_1*pixel) & (radius>=radius_2*pixel)
            exponential_1=zeros(N_filter,N_filter);
            for row=1:N_filter
                for column=1:N_filter
                    argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                    exponential_1(row,column)=exp(i*omega_0*argument);
                end
            end
            aerial_1=aerial_1+yita(p,q)* abs(  imfilter(double(m_binary),h.*exponential_1)  ).^2;
        end
    end
end 
z_binary=aerial_1>t_r;
sum6=sum(sum(abs(pz-double(z_binary))));   %Output pattern error of binary optimized mask
subplot(2,3,6);
imshow(z_binary);
xlabel(strcat('error=',num2str(sum6)));

%%%%%%Convergence of optimization algorithm%%%%%%
figure   
plot([1:count],convergence(1:count,:),'k');
title('Convergence of optimization algorithm');

%%%%%%Save all of the data%%%%%%
save Data_OPC_fse.mat;


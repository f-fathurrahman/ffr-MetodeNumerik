function [] = smo_PSM_mask(N, pz, N_filter, pixel, k, NA, lamda, order, sigma);
 
% clc;
% clear;
% %%%%%%The initialization of the parameter in the optimization%%%%%%
% N=80;   %Mask dimension
% N_filter=21;   %Amplitude impulse response dimension
% pixel=15;   %Pixel size (nano meter)
% k=0.29;
% NA=1.25;   %Numerical aperture
% lamda=193;   %Wavelength
% order=1;   %Order of Bessel function
% sigma=0.4;   %Partial coherence factor

midway=(N_filter+1)/2;   %middle of low pass filter

%%%%%%the amplitude impulse response of the partially coherent imaging system%%%%%%
h=zeros(N_filter,N_filter);
radius=0;
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
h=h/sum(sum(h));

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
 
% %%%%%%The desired output pattern%%%%%%
% pz=zeros(N,N);
% for ii=9:11
%     for jj=6:75
%         pz(ii,jj)=1;
%     end
% end
% for ii=15:17
%     for jj=17:32
%         pz(ii,jj)=-1;
%     end
% end
% for ii=18:20
%     for jj=17:24
%         pz(ii,jj)=-1;
%     end
% end
% for ii=21:23
%     for jj=9:24
%         pz(ii,jj)=-1;
%     end
% end
% for ii=27:29
%     for jj=9:24
%         pz(ii,jj)=1;
%     end
% end
% for ii=30:32
%     for jj=17:24
%         pz(ii,jj)=1;
%     end
% end
% for ii=33:35
%     for jj=17:32
%         pz(ii,jj)=1;
%     end
% end
% for ii=15:17
%     for jj=57:72
%         pz(ii,jj)=-1;
%     end
% end
% for ii=18:20
%     for jj=57:64
%         pz(ii,jj)=-1;
%     end
% end
% for ii=21:23
%     for jj=49:64
%         pz(ii,jj)=-1;
%     end
% end
% for ii=27:29
%     for jj=49:64
%         pz(ii,jj)=1;
%     end
% end
% for ii=30:32
%     for jj=57:64
%         pz(ii,jj)=1;
%     end
% end
% for ii=33:35
%     for jj=57:72
%         pz(ii,jj)=1;
%     end
% end
% for ii=39:41
%     for jj=6:75
%         pz(ii,jj)=-1;
%     end
% end
% for ii=45:47
%     for jj=17:32
%         pz(ii,jj)=1;
%     end
% end
% for ii=48:50
%     for jj=17:24
%         pz(ii,jj)=1;
%     end
% end
% for ii=51:53
%     for jj=9:24
%         pz(ii,jj)=1;
%     end
% end
% for ii=57:59
%     for jj=9:24
%         pz(ii,jj)=-1;
%     end
% end
% for ii=60:62
%     for jj=17:24
%         pz(ii,jj)=-1;
%     end
% end
% for ii=63:65
%     for jj=17:32
%         pz(ii,jj)=-1;
%     end
% end
% for ii=45:47
%     for jj=57:72
%         pz(ii,jj)=1;
%     end
% end
% for ii=48:50
%     for jj=57:64
%         pz(ii,jj)=1;
%     end
% end
% for ii=51:53
%     for jj=49:64
%         pz(ii,jj)=1;
%     end
% end 
% for ii=57:59
%     for jj=49:64
%         pz(ii,jj)=-1;
%     end
% end
% for ii=60:62
%     for jj=57:64
%         pz(ii,jj)=-1;
%     end
% end
% for ii=63:65
%     for jj=57:72
%         pz(ii,jj)=-1;
%     end
% end
% for ii=69:71
%     for jj=6:75
%         pz(ii,jj)=1;
%     end
% end
% figure
% imshow(pz,[-1,1]);
% axis on;

%%%%%%%Initial illumination pattern%%%%%%%
D=pixel*N;
D_C_1=lamda/2/sigma/NA;   %Coherence length
omega_0=pi/D;
N_coherence=18;   %Source dimension
midway_coherence=(N_coherence+1)/2;   %Middle point of illumination
radius_1=D/(2*D_C_1);   %Radius of annular illumination
yita=zeros(N_coherence,N_coherence);   %Illumination pattern
for row=1:N_coherence
    for column=1:N_coherence
        radius=pixel*sqrt( (row-midway_coherence)^2 + (column-midway_coherence)^2 );
        if (radius<=radius_1*pixel) 
            yita(row,column)=1;
        end
    end
end
figure
imshow(yita,[-1,1]);
axis on;
normalize=sum(sum(yita));   %Normalization factor

yita_initial=zeros(N_coherence,N_coherence);   %Initial illumination pattern
yita_initial=yita;


%%%%%%Initialization of optimization%%%%%%
threshold=0;
count=0;   %Iteration number
m=pz;   %Initial mask
direction_mask=zeros(N,N);   %Cost sensitivity function of mask
direction_source=zeros(N_coherence,N_coherence);   %Cost sensitivity function of source
flag_mask_pre=ones(N,N);   %Locations of the changable pixels on mask in the previous iteration
flag_source_pre=ones(N_coherence,N_coherence);   %Locations of the changable pixels on source in the previous iteration
flag_mask=zeros(N,N);   %Locations of the changable pixels on mask in the current iteration
flag_source=zeros(N_coherence,N_coherence);   %Locations of the changable pixels on source in the current iteration
error=10;   %Output pattern error in the current iteration
normalize=sum(sum(yita));

%%%%%%Calculate the output pattern error in the previous iteration%%%%%%
aerial=zeros(N,N);
for p=1:N_coherence
    for q=1:N_coherence    
        if (yita(p,q)>0)
            exponential=zeros(N_filter,N_filter); % has the same dimension as the filter h
            for row=1:N_filter
                for column=1:N_filter
                    argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                    exponential(row,column)=exp(i*omega_0*argument);
                end
            end
            aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
        end
    end
end
error_pre=sum(sum((abs(pz)-aerial).^2));   %Output pattern error in the previous iteration
disp(error_pre);
%%%%%%Source and mask optimization%%%%%%
while (error > threshold )   
   count=count+1;   %Update the iteration number
   disp(count);
   direction_mask=zeros(N,N);   %Reset the cost sensitivity function of mask
   direction_source=zeros(N_coherence,N_coherence);   %Reset the cost sensitivity function of source
   flag_mask=zeros(N,N);   %Reset the locations of the changable pixels on mask in the current iteration
   flag_source=zeros(N_coherence,N_coherence);   %Reset the locations of the changable pixels on source in the current iteration
   
   %%%%%%Calculate the aerial image%%%%%%
   aerial=zeros(N,N);
   normalize=sum(sum(yita));
   for p=1:N_coherence
       for q=1:N_coherence    
           if (yita(p,q)>0)
               exponential=zeros(N_filter,N_filter); % has the same dimension as the filter h
               for row=1:N_filter
                   for column=1:N_filter
                       argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                       exponential(row,column)=exp(i*omega_0*argument);
                   end
               end
               aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
           end
       end
   end
   
   %%%%%%Calculate the cost sensitivity function of mask%%%%%%
   for p=1:N_coherence
       for q=1:N_coherence   
           if (yita(p,q)>0)
               exponential=zeros(N_filter,N_filter);
               for row=1:N_filter
                   for column=1:N_filter
                       argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                       exponential(row,column)=exp(i*omega_0*argument);
                   end
               end
               direction_mask=direction_mask+(2)*yita(p,q)/normalize* real(  imfilter( (abs(pz)-aerial).*imfilter(double(m),h.*exponential) , conj(g.*exponential) )  +   imfilter( (abs(pz)-aerial).*imfilter(double(m),conj(h.*exponential)) , g.*exponential )    );%%%%%%%%%%%%get abs foe revision%%%
           end
       end
   end

   %%%%%%Calculate the locations of the changable pixels on mask%%%%%%
   for ii=2:N-1
       for jj=2:N-1
           if (m(ii,jj-1)+m(ii,jj+1)+m(ii-1,jj)+m(ii+1,jj)==0)|(m(ii,jj-1)+m(ii,jj+1)+m(ii-1,jj)+m(ii+1,jj)==4)
               flag_mask(ii,jj)=0;
           else
               flag_mask(ii,jj)=1;
           end
       end
   end
   
   %%%%%%If no pixel on mask or souce changed, then the program is terminated%%%%%%
   if(sum(sum(flag_mask~=flag_mask_pre)) == 0)
       break;
   else
       flag_mask_pre=flag_mask;
   end
   
   %%%%%%Search the changable pixels to reduce the output pattern error%%%%%%
   while (sum(sum(flag_mask))+sum(sum(flag_source))>0)   %Check whether there is any changable pixel on mask and source
       test_mask=abs(direction_mask.*flag_mask);   %Magnitude of the gradient for the changable pixels on mask
       max_index_mask=find(  test_mask >= max(max(test_mask))   );   %Find the changable pixel with maximum gradient magnitude on mask
       max_index_mask=max_index_mask(1,1);   %Choose a changable pixel with maximum gradient magnitude on mask
       max_row_mask=mod(max_index_mask,N);   %Determine the row of the changable pixel on mask
       if (max_row_mask==0)
           max_row_mask=N;
       end
       max_column_mask=floor(max_index_mask/N)+1;   %Determine the column of the changable pixel on mask
       
       %%%%%%Change the pixel on mask first%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (test_mask(max_row_mask,max_column_mask) >= 0)
           flag_mask(max_row_mask,max_column_mask)=0;   %Do not consider this pixel in the following iterations
           if ((m(max_row_mask,max_column_mask)<1)&(direction_mask(max_row_mask,max_column_mask)>0)) | ((m(max_row_mask,max_column_mask)>-1)&(direction_mask(max_row_mask,max_column_mask)<0))
               m(max_row_mask,max_column_mask)=m(max_row_mask,max_column_mask)+sign(direction_mask(max_row_mask,max_column_mask));   %Flip the pixel according to the cost sensitivity and current pixel value
               
               %%%%%%Check topological constraint on mask%%%%%%
               maxu=0;   %Flag of the topological constraint on mask
                         %If maxu=0, then the topological constraint on mask is satisfied
                         %If maxu=1, then the topological constraint on mask is not satisfied
               for ii=max(2,max_row_mask-1):min(N-1,max_row_mask+1)
                   if (maxu==1)
                       break;
                   end
                   for jj=max(2,max_column_mask-1):min(N-1,max_column_mask+1)
                       if (m(ii,jj)~=m(ii,jj-1)) & (m(ii,jj)~=m(ii,jj+1)) & (m(ii,jj)~=m(ii-1,jj)) & (m(ii,jj)~=m(ii+1,jj))
                           maxu=1;
                           break;
                       end
                   end
               end
               if (maxu==1)   %If the topological constraint on mask is not satisfied, then restore the original pixel values
                   m(max_row_mask,max_column_mask)=m(max_row_mask,max_column_mask)-sign(direction_mask(max_row_mask,max_column_mask));
                   continue;
               end
               
               %%%%%%Check the change of output pattern error%%%%%
               normalize=sum(sum(yita));
               aerial=zeros(N,N);
               for p=1:N_coherence
                   for q=1:N_coherence    
                       if (yita(p,q)>0)
                           exponential=zeros(N_filter,N_filter); % has the same dimension as the filter h
                           for row=1:N_filter
                               for column=1:N_filter
                                   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                                   exponential(row,column)=exp(i*omega_0*argument);
                               end
                           end
                           aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
                       end
                   end
               end
               error=sum(sum((abs(pz)-aerial).^2));   %Output pattern error in the current iteration
               
               if(error>=error_pre)   %If the error is increased, then restore the original pixel values
                   m(max_row_mask,max_column_mask)=m(max_row_mask,max_column_mask)-sign(direction_mask(max_row_mask,max_column_mask));
               else   %If the error is reduced, then keep the pixel update, and update the output patter error in the previous iteration
                   error_pre=error;
                   disp('change mask');
                   disp(error);
               end
           end         
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%The initial illumination%%%%%%
figure
subplot(2,3,1);
imshow(yita_initial,[-1,1]);
axis on;
title('Initial source');

%%%%%%The initial mask%%%%%%
subplot(2,3,2);
imshow(pz,[-1,1]);
axis on;
title('Initial mask');

%%%%%%The output pattern of initial mask and source%%%%%%
normalize=sum(sum(yita_initial));
aerial=zeros(N,N);
for p=1:N_coherence
    for q=1:N_coherence    
        if (yita_initial(p,q)>0)
           exponential=zeros(N_filter,N_filter);
           for row=1:N_filter
               for column=1:N_filter
                   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                   exponential(row,column)=exp(i*omega_0*argument);
               end
           end
           aerial=aerial+yita_initial(p,q)/ normalize* abs(  imfilter(double(pz),h.*exponential)  ).^2;
        end
    end
end
sum10=sum(sum((abs(pz)-aerial).^2));   %The output pattern error of initial mask and source

subplot(2,3,3);
imshow(aerial,[0,max(max(aerial))]);
colorbar('vertical');
xlabel(strcat('Error=',num2str(sum10)));
axis on;
title('Output');

%%%%%%Optimized illumination%%%%%%
subplot(2,3,4);
imshow(yita,[-1,1]);
axis on;
title('Optimized source');

%%%%%%The initial illumination%%%%%%
subplot(2,3,5);
imshow(m,[-1,1]);
axis on;
title('Optimized mask');

%%%%%%The output pattern of optimized mask and source%%%%%%
normalize=sum(sum(yita));
aerial=zeros(N,N);
for p=1:N_coherence
    for q=1:N_coherence    
        if (yita(p,q)>0)
           exponential=zeros(N_filter,N_filter);
           for row=1:N_filter
               for column=1:N_filter
                   argument=(p-midway_coherence)*(row-midway)*pixel+(q-midway_coherence)*(column-midway)*pixel;
                   exponential(row,column)=exp(i*omega_0*argument);
               end
           end
           aerial=aerial+yita(p,q)/ normalize* abs(  imfilter(double(m),h.*exponential)  ).^2;
        end
    end
end
sum6=sum(sum((abs(pz)-aerial).^2));   %The output pattern error of optimized mask and source

subplot(2,3,6);
imshow(aerial,[0,max(max(aerial))]);
colorbar('vertical');
xlabel(strcat('Error=',num2str(sum6)));
axis on;
title('Output');

%%%%%%Save all of the data%%%%%%
save Data_smo_PSM_mask.mat;



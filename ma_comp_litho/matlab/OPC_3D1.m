function [] = OPC_3D1(N, pz, N_filter, order);

% clc;
% clear;
% 
% N_filter=121;
% N=90;
% order=1;

N_dummy=N;
midway=(N_filter+1)/2;   %Middle of low pass filter
pixel=24.8;   %Pixel size
NA=0.68/4;
lamda=248;   %Wavelength

%%%%%%Optics low pass filter%%%%%%
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

% %%%%%%Desired output pattern%%%%%%
% pz=zeros(N,N);
% for ii=19:72
%     for jj=19:36
%         pz(ii,jj)=1;
%     end
% end
% 
% for ii=19:72
%     for jj=55:72
%         pz(ii,jj)=1;
%     end
% end

mask_initial=zeros(N,N);
for x=2:N-1
    for y=2:N-1
        if sum(sum(pz(x-1:x+1,y-1:y+1)))>0
            mask_initial(x,y)=1;
        end
    end
end
error_initial=sum(sum((pz-imfilter(pz,h).^2).^2));
disp(error_initial);

%%%%%%%Inialize other vectors%%%%%%
error_pre=0;
error=0;
gradient=zeros(N,N);
gradient_1=zeros(N,N);
m=zeros(N,N);
m=pz;
m_pre=zeros(N,N);
m_1_pre=zeros(8,8);
m_2_pre=zeros(8,8);

count=0;
count_small=0;
flag=0;
error_pre=sum(sum((pz-imfilter(m,h).^2).^2));
error=error_pre;

while ( sum(sum(m_pre~=m))>0 )
    
    m_pre=m;
    count=count+1;
    disp(count);
    gradient_1=(pz-imfilter(m,h).^2).*imfilter(m,h);
    gradient=-1*imfilter(gradient_1,g);

    for x=4.5:N-3.5
        for y=4.5:N-3.5
            m_1_pre=m(x-3.5:x+3.5,y-3.5:y+3.5);
            flag=0;
            if sum(sum(gradient(x-3.5:x+3.5,y-3.5:y+3.5)>0))>=8
                m(x-3.5:x+3.5,y-3.5:y+3.5)=0;   %Flip the pixel value 
                error=sum(sum((pz-imfilter(m,h).^2).^2));

                flag=check_OPC(m,N,8,3);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)
                    m(x-3.5:x+3.5,y-3.5:y+3.5)=m_1_pre;   %If the cost function is increased or the flag=1, then restore the pixel value                  
                else
                    error_pre=error;
                    disp('mask is changed');
                    disp(error);
                    disp(flag);
                end
            end

            m_1_pre=m(x-3.5:x+3.5,y-3.5:y+3.5);           
            flag=0;
            if sum(sum(gradient(x-3.5:x+3.5,y-3.5:y+3.5)<0))>=8
                m(x-3.5:x+3.5,y-3.5:y+3.5)=1;   %Flip the pixel value                              
                error=sum(sum((pz-imfilter(m,h).^2).^2));
                
                flag=check_OPC(m,N,8,3);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)
                    m(x-3.5:x+3.5,y-3.5:y+3.5)=m_1_pre;   %If the cost function is increased or the flag=1, then restore the pixel value                                   
                else
                    error_pre=error;
                    disp('mask is changed');
                    disp(error);
                    disp(flag);             
                end
            end          
        end
    end
end

figure
imshow(imfilter(m,h).^2,[0,1.1]);
title('Output pattern of the optimized binary mask');
axis on;
colorbar('vertical');

mask=zeros(N,N);
for x=2:N-1
    for y=2:N-1
        if sum(sum(m(x-1:x+1,y-1:y+1)))>0
            mask(x,y)=1;
        end
    end
end

figure
imshow(mask);
axis on;
title('Optimized binary mask');

save Data_OPC_3D1.mat;

 

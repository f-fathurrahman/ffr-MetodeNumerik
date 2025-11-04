function [] = OPC_3D2(N, pz, N_filter, order);

% clc;
% clear;
% 
% N_filter=139;
% N=95;
% order=1;

N_dummy=N;
midway=(N_filter+1)/2;   %Middle of low pass filter
pixel=14.5;   %Pixel size
NA=0.85/4;
lamda=193;   %Wavelength

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
% for ii=20:76
%     for jj=20:38
%         pz(ii,jj)=1;
%     end
% end
% 
% for ii=20:76
%     for jj=58:76
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

field_initial=pz;
for p=2:N-1
    for q=2:N-1
        if ((pz(p,q)==0) & (pz(p+1,q)==1)) | ((pz(p,q)==0) & (pz(p-1,q)==1))
            field_initial(p,q)=0.8i;
        end
    end
end
error_initial=sum(sum((pz-abs(imfilter(field_initial,h)).^2).^2));
disp(error_initial);

%%%%%%%Inialize other vectors%%%%%%
error_pre=0;
error=0;
m=zeros(N,N);
m=pz;
m_pre=zeros(N,N);
m_1_pre=zeros(12,12);
m_2_pre=zeros(12,12);
m_up=zeros(N,N);
m_down=zeros(N,N);
field=zeros(N,N);
count=0;
flag=0;
error_pre=sum(sum((pz-abs(imfilter(field_initial,h)).^2).^2));
error=error_pre;

while ( sum(sum(m_pre~=m))>0 )
    
    m_pre=m;
    count=count+1;
    disp(count);
    
    %%%%%%Shifting version of the mask pattern%%%%%%
    m_up=zeros(N,N);
    m_down=zeros(N,N);
    for x=1:N-1
        m_up(x,:)=m(x+1,:);
    end
    for x=2:N
        m_down(x,:)=m(x-1,:);
    end
    
    %%%%%%Caluculate the near field%%%%%%
    field=m;
    for p=2:N-1
        for q=2:N-1
            if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                field(p,q)=0.8i;
            end
        end
    end
    
    %%%%%%Calculate the cost sensitivity%%%%%%
    gradient=zeros(N,N);
    gradient_1=zeros(N,N);
    gradient_2_mid=zeros(N,N);
    gradient_2=zeros(N,N);
    gradient_3_mid=zeros(N,N);
    gradient_3=zeros(N,N);
    gradient_1=imfilter(   (pz-abs(imfilter(field,h)).^2).*imfilter(field,h),   g);
    gradient_2_mid=imfilter(   (pz-abs(imfilter(field,h)).^2).*imfilter(field,h),   g);
    for x=2:N
        gradient_2(x,:)=gradient_2_mid(x-1,:);
    end
    gradient_3_mid=imfilter(   (pz-abs(imfilter(field,h)).^2).*imfilter(field,h),   g);
    for x=1:N-1
        gradient_3(x,:)=gradient_3_mid(x+1,:);
    end
    gradient=-4*real( gradient_1.*(0.8i*m_up+0.8i*m_down+1)  +gradient_2.*(-0.8i).*(1-m_down)   +gradient_3.*(-0.8i).*(1-m_up) );

    %%%%%%Begin to scan the mask pattern to flip the right pixels%%%%%%
    for x=6.5:N-5.5
        for y=6.5:N-5.5
            m_1_pre=m(x-5.5:x+5.5,y-5.5:y+5.5);     
            flag=0;
            if sum(sum(gradient(x-5.5:x+5.5,y-5.5:y+5.5)>0))>=12
                m(x-5.5:x+5.5,y-5.5:y+5.5)=0;   %Flip the pixel value
                
                %%%%%%Calculate the near field%%%%%%
                field=m;
                for p=2:N-1
                    for q=2:N-1
                        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                            field(p,q)=0.8i;
                        end
                    end
                end
                error=sum(sum((pz-abs(imfilter(field,h)).^2).^2));

                flag=check_OPC(m,N,12,3);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)
                    m(x-5.5:x+5.5,y-5.5:y+5.5)=m_1_pre;   %If the cost function is increased or the flag=1, then restore the pixel value                  
                else
                    error_pre=error;
                    disp('mask is changed');
                    disp(error);
                    disp(flag);
                end
            end
   
            m_1_pre=m(x-5.5:x+5.5,y-5.5:y+5.5);  
            flag=0;
            if sum(sum(gradient(x-5.5:x+5.5,y-5.5:y+5.5)<0))>=12
                m(x-5.5:x+5.5,y-5.5:y+5.5)=1;   %Flip the pixel value
                
                %%%%%%Calculate the near field%%%%%%                
                field=m;
                for p=2:N-1
                    for q=2:N-1
                        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                            field(p,q)=0.8i;
                        end
                    end
                end
                error=sum(sum((pz-abs(imfilter(field,h)).^2).^2));

                flag=check_OPC(m,N,12,3);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)
                    m(x-5.5:x+5.5,y-5.5:y+5.5)=m_1_pre;   %If the cost function is increased or the flag=1, then restore the pixel value                    
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
imshow(m);
title('Optimized transmission area');

field=m;
for p=2:N-1
    for q=2:N-1
        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
            field(p,q)=0.8i;
        end
    end
end
figure
imshow(abs(field),[]);
title('Optimized near field');

figure
imshow(abs(imfilter(field,h)).^2,[0,2.4]);
axis on;
colorbar('vertical');
title('Output pattern of the optimized binary mask');

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

save Data_OPC_3D2.mat;

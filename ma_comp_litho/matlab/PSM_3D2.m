function [] = PSM_3D2(N, pz, N_filter, order);

% clc;
% clear;
% 
% N_filter=139;
% N=151;
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

 
%%%%%%Desired output pattern%%%%%%
% pz=zeros(N,N);
% for ii=48:104
%     for jj=48:66
%         pz(ii,jj)=1;
%     end
% end
% 
% for ii=48:104
%     for jj=86:104
%         pz(ii,jj)=-1;
%     end
% end




%%%%%%Initial mask%%%%%%
mask_initial_1=pz;
for x=5:N-4
    for y=5:N-4
        if sum(sum(    pz(x-4:x+4,y-4:y+4)<0   ))>0
            mask_initial_1(x,y)=-1;
        end
    end
end

mask_initial=mask_initial_1;
for x=2:N-1
    for y=2:N-1
        if sum(sum(    mask_initial_1(x-1:x+1,y-1:y+1)>0   ))>0
            mask_initial(x,y)=1;
        end
    end
end
figure
imshow(mask_initial,[-1,1]);
title('Initial PSM mask');

%%%%%%Initial near field%%%%%%
field_initial=pz;
for p=5:N-4
    for q=5:N-4
        if ((pz(p,q)==0) & (pz(p+4,q)==-1)) | ((pz(p,q)==0) & (pz(p-4,q)==-1))
            field_initial(p,q)=-0.3i;
        end
        if ((pz(p,q)==0) & (pz(p+1,q)==1)) | ((pz(p,q)==0) & (pz(p-1,q)==1))
            field_initial(p,q)=0.8i;
        end
    end
end
figure
imshow(real(field_initial)+imag(field_initial) ,[]);
title('Initial field');

%%%%%%Aerial image of initial mask%%%%%%
figure
imshow(abs(imfilter(field_initial,h)).^2,[]); 
axis on;
colorbar('vertical');
title('Output pattern of the initial mask');

%%%%%%Output error of initial mask%%%%%%
error_initial=sum(sum((abs(pz)-abs(imfilter(field_initial,h)).^2).^2));
disp(error_initial);

%%%%%%%Inialize other vectors%%%%%%
error_pre=0;
error=0;
m=zeros(N,N);
m=pz;
m_pre=zeros(N,N);
m_1_pre=zeros(16,16);
m_2_pre=zeros(16,16);
field_1_pre=zeros(24,24);
field_2_pre=zeros(24,24);

m_up=zeros(N,N);
m_up_4=zeros(N,N);
m_down=zeros(N,N);
m_down_4=zeros(N,N);
field=zeros(N,N);
count=0;
flag=0;
error_pre=sum(sum((abs(pz)-abs(imfilter(field_initial,h)).^2).^2));
error=error_pre;

while ( sum(sum(m_pre~=m))>0 )
    
    m_pre=m;
    count=count+1;
    disp(count);
    
    %%%%%%Shifting version of the mask pattern%%%%%%    
    m_up=zeros(N,N);
    m_up_4=zeros(N,N);
    m_down=zeros(N,N);
    m_down_4=zeros(N,N);
    for x=1:N-1
        m_up(x,:)=m(x+1,:);
    end
    for x=1:N-4
        m_up_4(x,:)=m(x+4,:);
    end
    for x=2:N
        m_down(x,:)=m(x-1,:);
    end
    for x=5:N
        m_down_4(x,:)=m(x-4,:);
    end
    
    %%%%%%Caluculate the near field%%%%%% 
    field=m;
    for p=5:N-4
        for q=5:N-4
            if ((m(p,q)==0) & (m(p+4,q)==-1)) | ((m(p,q)==0) & (m(p-4,q)==-1))
                field(p,q)=-0.3i;
            end
            if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                field(p,q)=0.8i;
            end
        end
    end
    
    %%%%%%Calculate the cost sensitivity%%%%%%
    gradient=zeros(N,N);
    gradient_1=zeros(N,N);
    gradient_1_trailer=zeros(N,N);
    gradient_2=zeros(N,N);
    gradient_2_trailer=zeros(N,N);
    gradient_3=zeros(N,N);
    gradient_3_trailer=zeros(N,N);
    gradient_4=zeros(N,N);
    gradient_4_trailer=zeros(N,N);
    gradient_5=zeros(N,N);
    gradient_5_trailer=zeros(N,N);
    
    gradient_1_trailer=0.3i.*m.*m_down_4.*(1-m_down_4) + 0.3i.*m.*m_up_4.*(1-m_up_4) + 0.8i.*m.*m_down.*(1-m_down) + 0.8i.*m.*m_up.*(1-m_up) + 1;
    gradient_2_trailer=(1-m_down).*(1+m_down).*(1+2*m).*(-0.4i);
    gradient_3_trailer=(1-m_up).*(1+m_up).*(1+2*m).*(-0.4i);
    gradient_4_trailer=(1-m_down_4).*(1+m_down_4).*(1-2*m).*(-0.15i);
    gradient_5_trailer=(1-m_up_4).*(1+m_up_4).*(1-2*m).*(-0.15i);
    gradient_1=imfilter(   (abs(pz)-abs(imfilter(field,h)).^2).*imfilter(field,h),   g);
    for x=2:N
        gradient_2(x,:)=gradient_1(x-1,:);
    end
    for x=1:N-1
        gradient_3(x,:)=gradient_1(x+1,:);
    end
    for x=5:N
        gradient_4(x,:)=gradient_1(x-4,:);
    end
    for x=1:N-4
        gradient_5(x,:)=gradient_1(x+4,:);
    end
    
    gradient=-4*real( gradient_1.*gradient_1_trailer + gradient_2.*gradient_2_trailer + gradient_3.*gradient_3_trailer + gradient_4.*gradient_4_trailer + gradient_5.*gradient_5_trailer );

    %%%%%%Begin to scan the mask pattern to flip the right pixels%%%%%%
    for x=16.5:N/2-7.5        
        for y=16.5:N-15.5                        
            m_1_pre=m(x-7.5:x+7.5,y-7.5:y+7.5);
            m_2_pre=m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5);
            field_1_pre=field(x-11.5:x+11.5,y-11.5:y+11.5);
            field_2_pre=field((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5);
           
            flag=0;
            if sum(sum(gradient(x-7.5:x+7.5,y-7.5:y+7.5)>0))>=16
                %%%%%%Flip the pixel value%%%%%%
                vic=min(min(m(x-7.5:x+7.5,y-7.5:y+7.5)));
                if vic==-1
                    m(x-7.5:x+7.5,y-7.5:y+7.5)=-1;  
                    m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5)=-1;
                else
                    m(x-7.5:x+7.5,y-7.5:y+7.5)=vic-1;
                    m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5)=vic-1;
                end
                %%%%%%Calculate the near field%%%%%%
                field(x-11.5:x+11.5,y-11.5:y+11.5)=m(x-11.5:x+11.5,y-11.5:y+11.5);
                for p=x-11.5:x+11.5
                    for q=y-11.5:y+11.5
                        if ((m(p,q)==0) & (m(p+4,q)==-1)) | ((m(p,q)==0) & (m(p-4,q)==-1))
                            field(p,q)=-0.3i;
                        end
                        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                            field(p,q)=0.8i;
                        end
                    end
                end
                field((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5)=m((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5);
                for p=(N+1-x)-11.5:(N+1-x)+11.5
                    for q=y-11.5:y+11.5
                        if ((m(p,q)==0) & (m(p+4,q)==-1)) | ((m(p,q)==0) & (m(p-4,q)==-1))
                            field(p,q)=-0.3i;
                        end
                        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                            field(p,q)=0.8i;
                        end
                    end
                end
                error=sum(sum((abs(pz)-abs(imfilter(field,h)).^2).^2));

                flag=check_PSM(m,N,16,1,4);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)   %If the cost function is increased or the flag=1, then restore the pixel value              
                    m(x-7.5:x+7.5,y-7.5:y+7.5)=m_1_pre;
                    m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5)=m_2_pre;
                    field(x-11.5:x+11.5,y-11.5:y+11.5)=field_1_pre;
                    field((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5)=field_2_pre;
                else
                    error_pre=error;
                    disp('mask is changed to the directio of -1');
                    disp(error);
                    disp(flag);
                end
            end
                     
            m_1_pre=m(x-7.5:x+7.5,y-7.5:y+7.5);
            m_2_pre=m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5);
            field_1_pre=field(x-11.5:x+11.5,y-11.5:y+11.5);
            field_2_pre=field((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5);
            
            flag=0;
            if sum(sum(gradient(x-7.5:x+7.5,y-7.5:y+7.5)<0))>=16
                %%%%%%Flip the pixel value%%%%%%
                vic=max(max(m(x-7.5:x+7.5,y-7.5:y+7.5)));
                if vic==1
                    m(x-7.5:x+7.5,y-7.5:y+7.5)=1;
                    m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5)=1;
                else
                    m(x-7.5:x+7.5,y-7.5:y+7.5)=vic+1;
                    m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5)=vic+1;
                end
                
                %%%%%%Calculate the near field%%%%%%
                field(x-11.5:x+11.5,y-11.5:y+11.5)=m(x-11.5:x+11.5,y-11.5:y+11.5);
                for p=x-11.5:x+11.5
                    for q=y-11.5:y+11.5
                        if ((m(p,q)==0) & (m(p+4,q)==-1)) | ((m(p,q)==0) & (m(p-4,q)==-1))
                            field(p,q)=-0.3i;
                        end
                        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                            field(p,q)=0.8i;
                        end
                    end
                end
                field((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5)=m((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5);
                for p=(N+1-x)-11.5:(N+1-x)+11.5
                    for q=y-11.5:y+11.5
                        if ((m(p,q)==0) & (m(p+4,q)==-1)) | ((m(p,q)==0) & (m(p-4,q)==-1))
                            field(p,q)=-0.3i;
                        end
                        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
                            field(p,q)=0.8i;
                        end
                    end
                end
                error=sum(sum((abs(pz)-abs(imfilter(field,h)).^2).^2));

                flag=check_PSM(m,N,16,1,4);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)   %If the cost function is increased or the flag=1, then restore the pixel value              
                    m(x-7.5:x+7.5,y-7.5:y+7.5)=m_1_pre;
                    m((N+1-x)-7.5:(N+1-x)+7.5,y-7.5:y+7.5)=m_2_pre;
                    field(x-11.5:x+11.5,y-11.5:y+11.5)=field_1_pre;
                    field((N+1-x)-11.5:(N+1-x)+11.5,y-11.5:y+11.5)=field_2_pre;
                else
                    error_pre=error;
                    disp('mask is changed to the direction of 1');
                    disp(error);
                    disp(flag);             
                end
            end              
        end
    end
end

figure
imshow(m,[-1,1]);
title('Optimized transmission area');

field=m;
for p=5:N-4
    for q=5:N-4
        if ((m(p,q)==0) & (m(p+4,q)==-1)) | ((m(p,q)==0) & (m(p-4,q)==-1))
            field(p,q)=-0.3i;
        end
        if ((m(p,q)==0) & (m(p+1,q)==1)) | ((m(p,q)==0) & (m(p-1,q)==1))
            field(p,q)=0.8i;
        end
    end
end
figure
imshow(real(field)+imag(field),[]);
title('Optimized near field');

figure
imshow(abs(imfilter(field,h)).^2,[-1,1]);
title('Output pattern of the optimized binary mask');

mask_1=m;
for x=5:N-4
    for y=5:N-4
        if sum(sum(    m(x-4:x+4,y-4:y+4)<0   ))>0
            mask_1(x,y)=-1;
        end
    end
end
mask=mask_1;
for x=2:N-1
    for y=2:N-1
        if sum(sum(    mask_1(x-1:x+1,y-1:y+1)>0   ))>0
            mask(x,y)=1;
        end
    end
end
figure
imshow(mask,[-1,1]);
axis on;
title('Optimized PSM mask');

save Data_PSM_3D2.mat;

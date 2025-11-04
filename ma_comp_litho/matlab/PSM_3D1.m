function [] = PSM_3D1(N, pz, N_filter, order);

% clc;
% clear;

% N_filter=109;
% N=80;
% order=1;

N_dummy=N;
midway=(N_filter+1)/2;   %Middle of low pass filter
pixel=55.8/2;   %Pixel size
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
% for ii=17:64
%     for jj=17:32
%         pz(ii,jj)=1;
%     end
% end
% 
% for ii=17:64
%     for jj=49:64
%         pz(ii,jj)=-1;
%     end
% end

%%%%%%Initial mask%%%%%%
mask_initial_1=pz;
for x=3:N-2
    for y=3:N-2
        if sum(sum(    pz(x-2:x+2,y-2:y+2)<0   ))>0
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
axis on;
title('Initial mask');

%%%%%%Initial near field%%%%%%
field_initial=pz;
for p=3:N-2
    for q=3:N-2
        if ((pz(p,q)==0) & (pz(p+2,q)==-1)) | ((pz(p,q)==0) & (pz(p-2,q)==-1))
            field_initial(p,q)=-0.52i;
        end
    end
end
figure
imshow(real(field_initial)+imag(field_initial),[]);
axis on;
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
m_1_pre=zeros(9,9);
m_2_pre=zeros(9,9);
field_1_pre=zeros(13,13);
field_2_pre=zeros(13,13);
m_up=zeros(N,N);
m_down=zeros(N,N);
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
    m_down=zeros(N,N);
    for x=1:N-2
        m_up(x,:)=m(x+2,:);
    end
    for x=3:N
        m_down(x,:)=m(x-2,:);
    end
    
    %%%%%%Caluculate the near field%%%%%%    
    field=m;
    for p=3:N-2
        for q=3:N-2
            if ((m(p,q)==0) & (m(p+2,q)==-1)) | ((m(p,q)==0) & (m(p-2,q)==-1))
                field(p,q)=-0.52i;
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
    
    gradient_1_trailer=0.52i.*m.*m_up.*(1-m_up) + 0.52i.*m.*m_down.*(1-m_down) + 1;
    gradient_2_trailer=-0.26i.*(1-m_down).*(1+m_down).*(1-2.*m);
    gradient_3_trailer=-0.26i.*(1-m_up).*(1+m_up).*(1-2.*m);
    gradient_1=imfilter(   (abs(pz)-abs(imfilter(field,h)).^2).*imfilter(field,h),   g);
    for x=3:N
        gradient_2(x,:)=gradient_1(x-2,:);
    end
    for x=1:N-2
        gradient_3(x,:)=gradient_1(x+2,:);
    end
    
    gradient=-4*real( gradient_1.*gradient_1_trailer + gradient_2.*gradient_2_trailer + gradient_3.*gradient_3_trailer );

    %%%%%%Begin to scan the mask pattern to flip the right pixels%%%%%%
    for x=9:N/2-8
        for y=9:N-8                        
            m_1_pre=m(x-4:x+4,y-4:y+4);
            m_2_pre=m((N+1-x)-4:(N+1-x)+4,y-4:y+4);
            field_1_pre=field(x-6:x+6,y-6:y+6);
            field_2_pre=field((N+1-x)-6:(N+1-x)+6,y-6:y+6);
            
            flag=0;
            if sum(sum(gradient(x-4:x+4,y-4:y+4)>0))>=9
                %%%%%%Flip the pixel value%%%%%%
                vic=min(min(m(x-4:x+4,y-4:y+4)));
                if vic==-1
                    m(x-4:x+4,y-4:y+4)=-1;  
                    m((N+1-x)-4:(N+1-x)+4,y-4:y+4)=-1;
                else
                    m(x-4:x+4,y-4:y+4)=vic-1;
                    m((N+1-x)-4:(N+1-x)+4,y-4:y+4)=vic-1;
                end
                %%%%%%Calculate the near field%%%%%%                
                field(x-6:x+6,y-6:y+6)=m(x-6:x+6,y-6:y+6);
                for p=x-6:x+6
                    for q=y-6:y+6
                        if ((m(p,q)==0) & (m(p+2,q)==-1)) | ((m(p,q)==0) & (m(p-2,q)==-1))
                            field(p,q)=-0.52i;
                        end
                    end
                end
                field((N+1-x)-6:(N+1-x)+6,y-6:y+6)=m((N+1-x)-6:(N+1-x)+6,y-6:y+6);
                for p=(N+1-x)-6:(N+1-x)+6
                    for q=y-6:y+6
                        if ((m(p,q)==0) & (m(p+2,q)==-1)) | ((m(p,q)==0) & (m(p-2,q)==-1))
                            field(p,q)=-0.52i;
                        end
                    end
                end
                error=sum(sum((abs(pz)-abs(imfilter(field,h)).^2).^2));

                flag=check_PSM(m,N,9,1,2);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)   %If the cost function is increased or the flag=1, then restore the pixel value                  
                    m(x-4:x+4,y-4:y+4)=m_1_pre;
                    m((N+1-x)-4:(N+1-x)+4,y-4:y+4)=m_2_pre;
                    field(x-6:x+6,y-6:y+6)=field_1_pre;
                    field((N+1-x)-6:(N+1-x)+6,y-6:y+6)=field_2_pre;
                else
                    error_pre=error;
                    disp('mask is changed to the directio of -1');
                    disp(error);
                    disp(flag);
                end
            end
            

            
            
            m_1_pre=m(x-4:x+4,y-4:y+4);
            m_2_pre=m((N+1-x)-4:(N+1-x)+4,y-4:y+4);
            field_1_pre=field(x-6:x+6,y-6:y+6);
            field_2_pre=field((N+1-x)-6:(N+1-x)+6,y-6:y+6);
            
            flag=0;
            if sum(sum(gradient(x-4:x+4,y-4:y+4)<0))>=9
                %%%%%%Flip the pixel value%%%%%%
                vic=max(max(m(x-4:x+4,y-4:y+4)));
                if vic==1
                    m(x-4:x+4,y-4:y+4)=1;
                    m((N+1-x)-4:(N+1-x)+4,y-4:y+4)=1;
                else
                    m(x-4:x+4,y-4:y+4)=vic+1;
                    m((N+1-x)-4:(N+1-x)+4,y-4:y+4)=vic+1;
                end
                %%%%%%Calculate the near field%%%%%%  
                field(x-6:x+6,y-6:y+6)=m(x-6:x+6,y-6:y+6);
                for p=x-6:x+6
                    for q=y-6:y+6
                        if ((m(p,q)==0) & (m(p+2,q)==-1)) | ((m(p,q)==0) & (m(p-2,q)==-1))
                            field(p,q)=-0.52i;
                        end
                    end
                end
                field((N+1-x)-6:(N+1-x)+6,y-6:y+6)=m((N+1-x)-6:(N+1-x)+6,y-6:y+6);
                for p=(N+1-x)-6:(N+1-x)+6
                    for q=y-6:y+6
                        if ((m(p,q)==0) & (m(p+2,q)==-1)) | ((m(p,q)==0) & (m(p-2,q)==-1))
                            field(p,q)=-0.52i;
                        end
                    end
                end
                error=sum(sum((abs(pz)-abs(imfilter(field,h)).^2).^2));

                flag=check_PSM(m,N,9,1,2);   %Check the topological constraint

                if (error>=error_pre) | (flag==1)   %If the cost function is increased or the flag=1, then restore the pixel value                  
                    m(x-4:x+4,y-4:y+4)=m_1_pre;
                    m((N+1-x)-4:(N+1-x)+4,y-4:y+4)=m_2_pre;
                    field(x-6:x+6,y-6:y+6)=field_1_pre;
                    field((N+1-x)-6:(N+1-x)+6,y-6:y+6)=field_2_pre;
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
for p=3:N-2
    for q=3:N-2
        if ((m(p,q)==0) & (m(p+2,q)==-1)) | ((m(p,q)==0) & (m(p-2,q)==-1))
             field(p,q)=-0.52i;
        end
    end
end
figure
imshow(real(field)+imag(field),[]);
title('Optimized near field');

figure
imshow(abs(imfilter(field,h)).^2,[]);
axis on;
colorbar('vertical');
title('Output pattern of the optimized PSM mask');

mask_1=m;
for x=3:N-2
    for y=3:N-2
        if sum(sum(    m(x-2:x+2,y-2:y+2)<0   ))>0
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
title('Optimized mask');

save Data_PSM_3D1.mat;

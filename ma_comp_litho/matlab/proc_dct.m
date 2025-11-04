%%%%%%Post-processing based on 2D-DCT%%%%%%
function [m_trinary_new] = proc_dct(N_mask, pz, m, t_r, t_r_real, t_m, TCC, threshold);

error=zeros(101,1);  %The output pattern error versus the number of reserved DCT low frequency components

figure
%%%%%%Gray optimized mask%%%%%%
subplot(1,2,1);
imshow(m,[-1,1]);
axis on;
title('Gray mask');
%%%%%%2D-DCT of gray optimized mask%%%%%%
B=dct2(m);   %2D-DCT of gray optimized mask
subplot(1,2,2);
imshow(B,[-1,1]);
axis on;
title('2D-DCT');

B_new=B;   %2D-DCT transform after the cutting off the low frequency components
for ii=102:-1:2
    disp(ii);
    for p=1:N_mask
        for q=1:N_mask
            if (p+q==ii)
                B_new(p,q)=0;
            end
        end
    end
    m_new= idct2(B_new);   %Gray mask after the cutting off the low frequency components
    m_trinary_p=m_new>t_m;
    m_trinary_n=-1*(m_new<(-1*t_m));
    m_trinary_new=m_trinary_p+m_trinary_n;   %Trinary mask after the cutting off the low frequency components
    
    %%%%Calculate the output pattern error corresponding to the trinary%%%%
    %%%%mask after the cutting off the low frequency components%%%%%%%%%%%%
    aerial=zeros(N_mask,N_mask);
    aerial_fre=zeros(N_mask,N_mask);
    m_trinary_new_fre=(fftshift(fft2(m_trinary_new)));
    for x=1:N_mask^2
        for y=1:N_mask^2
            index_1=mod(x-1,N_mask)+1;
            index_2=floor((x-1)/N_mask)+1;
            index_3=mod(y-1,N_mask)+1;
            index_4=floor((y-1)/N_mask)+1;
            aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)=aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)+TCC(x,y)*(m_trinary_new_fre(index_1,index_2))*conj(m_trinary_new_fre(index_3,index_4));
        end
    end
    aerial=abs(ifft2(aerial_fre))/((N_mask)^2);
    z_trinary=aerial>t_r_real;   %Binary output pattern after cutting off the low frequency components
    error(103-ii,1)=sum(sum(abs(abs(pz)-z_trinary)));   %Output pattern error after cutting off the low frequency components
    disp(error(103-ii,1));
end

%%%%%%Display the relationship between the number of maintained DCT low%%%%
%%%%%%frequency components and the output pattern errors%%%%%%%%%%%%%%%%%%%

vic=ones(N_mask,N_mask);   %The location where DCT coefficients are maintained
curve=zeros(N_mask^2,1);   %The output pattern error with respect to the number of maintained DCT coefficients
for ii=102:-1:2
    for p=1:N_mask
        for q=1:N_mask
            if (p+q==ii)
                vic(p,q)=0;
            end
        end
    end
    curve(sum(sum(vic))+1)=error(103-ii,1);   %Interpolation of the curve
end
for ii=2:N_mask^2
    if (curve(ii)==0)
        curve(ii)=curve(ii-1);
    end
end
figure
plot([0:N_mask^2-1],curve);
axis([0 2600 0 500]);

%%%%%%The optimal post-processing chosen according to the relationship%%%%%%
B_final= B;   %2D-DCT transform after the optimal post-processing
for p=1:N_mask
    for q=1:N_mask
        if (p+q>=threshold)   %Maintain 136 low frequency components 
            B_final(p,q)=0;
        end
    end
end
m_new=idct2(B_final);   %Gray mask after the optimal post-processing
m_trinary_p=m_new>t_m;
m_trinary_n=-1*(m_new<(-1*t_m));
m_trinary_new=m_trinary_p+m_trinary_n;   %Trinary mask after the optimal post-processing

%%%%%%Display the trinary mask after the optimal post-processing%%%%%%
figure
imshow(m_trinary_new,[-1,1]);
title('Optimized mask');
axis on;
  
%%%%Calculate the output pattern error corresponding to the trinary%%%%
%%%%mask after the optimal post-processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aerial=zeros(N_mask,N_mask);
aerial_fre=zeros(N_mask,N_mask);
m_trinary_new_fre=(fftshift(fft2(m_trinary_new)));
for x=1:N_mask^2
    for y=1:N_mask^2
        index_1=mod(x-1,N_mask)+1;
        index_2=floor((x-1)/N_mask)+1;
        index_3=mod(y-1,N_mask)+1;
        index_4=floor((y-1)/N_mask)+1;
        aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)=aerial_fre(mod(index_1-index_3,N_mask)+1,mod(index_2-index_4,N_mask)+1)+TCC(x,y)*(m_trinary_new_fre(index_1,index_2))*conj(m_trinary_new_fre(index_3,index_4));
    end
    disp(x);
end
aerial=abs(ifft2(aerial_fre))/((N_mask)^2);
z_trinary=aerial>t_r_real;   %Binary output pattern after the optimal post-processing
error_final=sum(sum(abs(abs(pz)-z_trinary)));   %Output pattern error after the optimal post-processing

%%%%%%Display the output pattern after the optimal post-processing%%%%%%
figure
imshow(z_trinary,[-1,1]);
xlabel(strcat('Error=',num2str(error_final)));
title('Output pattern');
axis on;
 

 

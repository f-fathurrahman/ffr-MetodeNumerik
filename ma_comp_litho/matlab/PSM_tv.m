function [] = PSM_tv(N, N_filter, k, pz, r, s, a, t_r, t_m, gamma_D, gamma_TV, epsilon, maxloop);

% clc;
% clear;
%%%%%%The initialization of the parameter in the optimization%%%%%%
% N=50;   %Mask dimension
% s=1;   %Step size
% a=90;
% t_r=0.5;
% t_m=0.5;   %Global threshold of the mask
% gamma_D=0;   %Weight of the discretization penalty
% gamma_TV=0;   %Weight of the total variation penalty
% epsilon=10;   %Tolerable output pattern error
% maxloop=160;   %Maximum iteration number

d=zeros(N,N);   %Gradient of the cost function
d_D=zeros(N,N);   %Gradient of the discretization penalty
d_TV=zeros(N,N);   %Gradient of the total variation penalty
convergence=zeros(maxloop,1);   %Output pattern error in each iteration
count=0;   %Index of iteration number
sum6=100;   %Output pattern error corresponding to the optimized binary mask
sum8=100;   %Output pattern error corresponding to the optimized real-valued mask

%%%%%%the amplitude impulse response of the coherent imaging system%%%%%%
h=fspecial('gaussian',N_filter,k);

%%%%%%The rotation of the amplitude impulse response%%%%%%
for i=1:N_filter
    for j=1:N_filter
        h1((i-1)*N_filter+j)=h(i,j);
    end
end
for i=1:N_filter
    for j=1:N_filter
        g(i,j)=h1((N_filter-i)*N_filter+(N_filter+1-j));
    end
end

%%%%%%The desired output pattern%%%%%%
% pz=zeros(N,N);
% for i=16:35
%     for j=12:22
%         pz(i,j)=1;
%     end
% end
% for i=16:35
%     for j=28:38
%         pz(i,j)=1;
%     end
% end
figure
imshow(pz,[-1,1]);
axis on;

%%%%%%The initialization of \theta, where r=\theta%%%%%%
% r=ones(N,N)*pi/2;
% for i=16:35
%     for j=12:22
%         r(i,j)=0;
%     end
% end
% for i=16:35
%     for j=28:38
%         r(i,j)=pi;
%     end
% end

%%%%%%PSM optimization in coherent imaging system%%%%%%
m=zeros(N,N);   %Mask pattern
while (sum6>epsilon) & (count<maxloop)
    count=count+1; 
    r=r-s*d;   %Update
    m=cos(r);   %Update
    viccinone=m>t_m;
    viccintwo=m<-1*t_m;
    viccin=viccinone+(-1)*viccintwo;   %Trinary mask
    viccout=imfilter(viccin,h);
    viccbin=abs(viccout)>t_r;   %Output pattern of trinary mask
    sum6=sum(sum(abs(pz-viccbin)));  
    convergence(count,1)=sum6; 

    mid1=imfilter(m,h);
   
    compone=mid1>=0;
    comptwo=mid1<0;
    comp=compone+comptwo*-1;
   
    z=1./ (  1+exp(-1*a*abs(mid1)+a*t_r)  ); 
    mid3=( pz-z ).*z.*(1-z).*comp;      
    mid5=imfilter(mid3,g);

    %%%%%%Gradient of discretization penaly%%%%%%  
    d_D=2*m.*(-1*sin(r));
   
    %%%%%%Gradient of total variation penaly%%%%%%
    f=abs(abs(m)-pz);
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

    d_TV=(f11+f22).*sign(m-pz)*(-1).*sin(r);
   
    %%%%%%Gradient of overall cost function%%%%%%
    d=2*a*mid5.*sin(r) +gamma_D*d_D +gamma_TV*d_TV;
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Display%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Output pattern of desired pattern%%%%%%
output_pz=abs(imfilter(pz,h))>t_r;
figure   %Output pattern of desired pattern
imshow(output_pz,[-1,1]);
axis on;
sum10=sum(sum(abs(pz-output_pz)));   %Output pattern error of desired pattern

%%%%%%Real-valued optimized mask and its output pattern%%%%%%
figure   %Real-valued optimized mask
imshow(m,[-1,1]);
axis on;

output_m=abs(imfilter(m,h))>t_r;
figure   %Output pattern of real-valued optimized mask
imshow(output_m,[-1,1]);
axis on;
sum8=sum(sum(abs(pz-output_m)));   %Output pattern error of real-valued optimized mask

%%%%%%Binary optimized mask and its output pattern%%%%%%
m2one=m>t_m;
m2two=m<-1*t_m;
m_b=m2one+(-1)*m2two;
figure   %Binary optimized mask
imshow(m_b,[-1,1]);
axis on;

output_m_b=abs(imfilter(m_b,h))>t_r;
figure   %Output pattern of binary optimized mask 
imshow(output_m_b,[-1,1]);
axis on;
sum6=sum(sum(abs(pz-output_m_b)));   %Output pattern error of binary optimized mask 

figure   %Convergence of optimization algorithm
plot([1:count],convergence(1:count,:),'k');

%%%%%%Save all of the data%%%%%%
save Data_PSM_tv.mat;










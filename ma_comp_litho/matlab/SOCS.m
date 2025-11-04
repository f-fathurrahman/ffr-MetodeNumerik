%%%%%%Calculate the transmission cross coefficients%%%%%%
function [TCC] = SOCS(N, pixel, k, NA, lamda, midway, sigma, order);

J=zeros(N,N);   %Effective source
P=zeros(N,N);   %Pupil function = fourier transform of lens shape
TCC=zeros(N^2,N^2);   %Transmission cross-coefficient

%%%%%%Pupil function = fourier transform of lens shape%%%%%%
P=zeros(N,N);
radius_filter=NA/lamda*pixel*N;
for row=1:N
    for column=1:N
        if (radius_filter>=sqrt( (row-midway)^2 + (column-midway)^2 ) );
            P(row,column)=1/pi/radius_filter^2;
        end
    end
end

%%%%%%Amplitude impulse response%%%%%%
h=zeros(N,N);
h=fftshift(ifft2(ifftshift(P)));

%%%%%%Effective source%%%%%%
J=zeros(N,N);
radius_illumination=radius_filter*sigma;
for row=1:N
    for column=1:N
        if (radius_illumination>=sqrt( (row-midway)^2 + (column-midway)^2 ) );
            J(row,column)=1/pi/radius_illumination^2;
        end
    end
end

%%%%%%Calculate transmission cross coefficients%%%%%%
P_1=P;
P_2=P;
for x=1:N^2
    for y=1:N^2
        index_1=mod(x-1,N)+1-midway;
        index_2=floor((x-1)/N)+1-midway;
        index_3=mod(y-1,N)+1-midway;
        index_4=floor((y-1)/N)+1-midway;
        if (index_1>0)
            P_1=[P(index_1+1:N,1:N);zeros(index_1,N)];
        end
        if (index_1<0)
            P_1=[zeros(abs(index_1),N);P(1:N+index_1,1:N)];
        end
        if (index_2>0)
            P_1=[P_1(1:N,index_2+1:N),zeros(N,index_2)];
        end
        if (index_2<0)
            P_1=[zeros(N,abs(index_2)),P_1(1:N,1:N+index_2)];
        end


        if (index_3>0)
            P_2=[P(index_3+1:N,1:N);zeros(index_3,N)];
        end
        if (index_3<0)
            P_2=[zeros(abs(index_3),N);P(1:N+index_3,1:N)];
        end
        if (index_4>0)
            P_2=[P_2(1:N,index_4+1:N),zeros(N,index_4)];
        end
        if (index_4<0)
            P_2=[zeros(N,abs(index_4)),P_2(1:N,1:N+index_4)];
        end
        TCC(x,y)=sum(sum(J.*P_1.*P_2));
    end
    disp(x);
end




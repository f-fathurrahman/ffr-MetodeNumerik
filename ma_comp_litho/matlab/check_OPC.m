function flag = check_OPC(m_dummy,N_dummy,singular1,singular2);

symbol_1=m_dummy;
d1=(singular1+1)/2;
d2=(singular1-1)/2;
for x_dummy=d1:N_dummy-d2
    for y_dummy=d1:N_dummy-d2
        if sum(sum(m_dummy(x_dummy-d2:x_dummy+d2,y_dummy-d2:y_dummy+d2)))==(singular1)^2
            symbol_1(x_dummy-d2:x_dummy+d2,y_dummy-d2:y_dummy+d2)=0;
        end
    end
end
c1=(singular2+1)/2;
c2=(singular2-1)/2;
symbol_2=1-m_dummy;
for x_dummy=c1:N_dummy-c2
    for y_dummy=c1:N_dummy-c2
        if sum(sum(m_dummy(x_dummy-c2:x_dummy+c2,y_dummy-c2:y_dummy+c2)))==0
            symbol_2(x_dummy-c2:x_dummy+c2,y_dummy-c2:y_dummy+c2)=0;
        end
    end
end
if sum(sum(symbol_1))+sum(sum(symbol_2))>0
    flag=1;
else
    flag=0;
end

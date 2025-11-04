function flag = check_PSM(m_dummy,N_dummy,singular12,boundary_clear,boundary_shift);

symbol_1=double(m_dummy>=1);
d1=(singular12+1)/2;
d2=(singular12-1)/2;
for x_dummy=d1:N_dummy-d2
    for y_dummy=d1:N_dummy-d2
        if sum(sum(m_dummy(x_dummy-d2:x_dummy+d2,y_dummy-d2:y_dummy+d2)))==(singular12)^2
            symbol_1(x_dummy-d2:x_dummy+d2,y_dummy-d2:y_dummy+d2)=0;
        end
    end
end

symbol_3=double(m_dummy<=-1);
for x_dummy=d1:N_dummy-d2
    for y_dummy=d1:N_dummy-d2
        if sum(sum(m_dummy(x_dummy-d2:x_dummy+d2,y_dummy-d2:y_dummy+d2)))==-1*(singular12)^2
            symbol_3(x_dummy-d2:x_dummy+d2,y_dummy-d2:y_dummy+d2)=0;
        end
    end
end

singular3=max(boundary_clear,boundary_shift);
c1=singular3+1;
c2=singular3;
symbol_2=double(m_dummy==0);
for x_dummy=c1:N_dummy-c2
    for y_dummy=c1:N_dummy-c2
        if sum(sum(    abs(m_dummy(x_dummy-c2:x_dummy+c2,y_dummy-c2:y_dummy+c2))     ))==0
            symbol_2(x_dummy-c2:x_dummy+c2,y_dummy-c2:y_dummy+c2)=0;
        end
    end
end

singular4=boundary_clear+boundary_shift+1;
test=double(m_dummy<=-1);
maxu=double(m_dummy<=-1);
for x_dummy=singular4+1:N_dummy-singular4
    for y_dummy=singular4+1:N_dummy-singular4
        if sum(sum(    maxu(x_dummy-singular4:x_dummy+singular4,y_dummy-singular4:y_dummy+singular4)   ))>0
            test(x_dummy,y_dummy)=1;
        end
    end
end
symbol_4=(test+double(m_dummy>=1))>1;

if sum(sum(symbol_1))+sum(sum(symbol_2))+sum(sum(symbol_3))+sum(sum(symbol_4))>0
    flag=1;
else
    flag=0;
end


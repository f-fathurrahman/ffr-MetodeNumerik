function [dw,P,P_inc,r_M,ErrorM]  ...
                    = InitialEstimatedw(eta_M,data,input)
global nDIM 

NS = input.NS;  
dw = cell(nDIM,NS); P = cell(nDIM,NS);   P_inc = cell(nDIM,NS); 
r_M = cell(1,NS);

Norm_M = 0; 
dummyq = cell(1,nDIM);
for q = 1 : NS
    r_M{q}    = data(:,q);
    g_dw      = AdjDOPMdw(r_M{q},input); 
    A0        = norm(g_dw{1}(:))^2 + norm(g_dw{2}(:))^2;
    dummyq{1} = g_dw{1};
    dummyq{2} = g_dw{2};
    GMp       = DOPMdw(dummyq,input);
    B0        = norm(GMp(:))^2; 
    alpha     = real(A0/B0);
    r_M{q}    = r_M{q} - alpha * GMp(:);
    Norm_M    = Norm_M + norm(r_M{q})^2; 
    dw{1,q}   = alpha * g_dw{1};
    dw{2,q}   = alpha * g_dw{2}; 
    dummyq{1} = dw{1,q};
    dummyq{2} = dw{2,q};
         Pinc = IncPressureWave(q,input);
   P_inc{1,q} = Pinc{1};
   P_inc{2,q} = Pinc{2};
   
   if (input.Kirchhoff == 0)  || (input.Kirchhoff == 1)
      [Kdw]  = KOPdw(dummyq,input); 
   elseif input.Kirchhoff == 2
      Kdw{1} = 0 * dw{1,q};                      
      Kdw{2} = 0 * dw{2,q};  
    end  
      P{1,q} = P_inc{1,q} + Kdw{1};
      P{2,q} = P_inc{2,q} + Kdw{2};
end %qloop  

ErrorM = sqrt(eta_M * Norm_M);
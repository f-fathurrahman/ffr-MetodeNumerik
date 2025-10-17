function [wZv,Zv,r_M,ErrorM]  ...
                    = InitialEstimatewZv(eta_M,Zvinc,dG_R,data,input)
global nDIM 
NS = input.NS;  
wZv = cell(nDIM,NS); Zv = cell(nDIM,NS);  r_M = cell(1,NS);

Norm_M = 0;
dummyq = cell(1,nDIM); 
for q = 1 : NS
    r_M{q}    = data(:,q);  
    [gZv]     = AdjDOPMwZv(dG_R,r_M{q},input);                 
    A0        = norm(gZv{1}(:))^2 + norm(gZv{2}(:))^2 ;
    GMp       = DOPMwZv(dG_R,gZv,input);            
    B0        = norm(GMp(:))^2; 
    alpha     = real(A0/B0);
    r_M{q}    = r_M{q} - alpha * GMp(:);
    Norm_M    = Norm_M + norm(r_M{q})^2;     
    wZv{1,q}  = alpha * gZv{1};
    wZv{2,q}  = alpha * gZv{2};
    dummyq{1} = wZv{1,q}; 
    dummyq{2} = wZv{2,q}; 
    [KZv]     = KOPwZv(dummyq,input);
    Zv{1,q}   = Zvinc{1,q} + KZv{1};
    Zv{2,q}   = Zvinc{2,q} + KZv{2};
end %qloop 

ErrorM = sqrt(eta_M * Norm_M);
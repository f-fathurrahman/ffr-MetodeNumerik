function [wZv] = ITERCGMDwZv(Zvinc,dG_R,data,input) 
global nDIM it
   NS = input.NS;
  AN1 = cell(1,NS);        
 wZvq = cell(nDIM);      vZvq = cell(nDIM);      Zvq = cell(nDIM);    
g1Zvq = cell(nDIM);     rDZvq = cell(nDIM);      
g1Zv  = cell(nDIM,NS);   vZv  = cell(nDIM,NS);

itmax = 256;   
eta_M = 1 / norm(data(:))^2;     % normalization factor eta_M
 
% (1) Initial estimate of contrast sources and  mass-density contrast
   [wZv,Zv,r_M,ErrorM] = InitialEstimatewZv(eta_M,Zvinc,dG_R,data,input);
   it = 0;  ErrorD  = 1;   
   disp(['Iteration = ', num2str(it)]);
   [CHI_rho,rZv,eta_D,ErrorD] ...
                   = UpdateContrastwZv(wZv,Zv,Zvinc,ErrorM,ErrorD,input);       
   saveCHI(CHI_rho); 
   disp('--------------------------------------------------------------');
 
% (2) Iterative updating of contrast sources and mass-density contrast  
   it = 1;  
   while (it <= itmax)    
     Norm_M = 0;
     for q = 1 : NS  
       AN1q = AN1{q}; 
       for n = 1 : nDIM 
        wZvq{n} =  wZv{n,q};   vZvq{n} = vZv{n,q};    Zvq{n} = Zv{n,q};
       g1Zvq{n} = g1Zv{n,q};    rMq    = r_M{q};    rDZvq{n} = rZv{n,q};
      end 
      [wZvq,vZvq,Zvq,g1Zvq,AN1q,rMq] = InvITERwZv(CHI_rho,dG_R,wZvq, ...
                         vZvq,Zvq,g1Zvq,AN1q,rMq,rDZvq,eta_M,eta_D,input);           
      AN1{q} = AN1q;   
      for n = 1 : nDIM 
        wZv{n,q} =  wZvq{n};  vZv{n,q} = vZvq{n};  Zv{n,q} = Zvq{n};
       g1Zv{n,q} = g1Zvq{n};    r_M{q} = rMq;
      end 
      Norm_M = Norm_M + norm(r_M{q})^2; 
    end % qloop   
    ErrorM = eta_M * Norm_M;  
    disp(['Iteration = ', num2str(it)]); 
    [CHI_rho,rZv,eta_D,ErrorD] ...
                   = UpdateContrastwZv(wZv,Zv,Zvinc,ErrorM,ErrorD,input);
    saveCHI(CHI_rho);
    it = it+1;
   end % while 
disp('--------------------------------------------------------------');
Error  = sqrt(ErrorM + ErrorD);  
disp(['Number of iterations = ' num2str(it)]);
disp(['Total Error in M and D = ' num2str(Error)]);
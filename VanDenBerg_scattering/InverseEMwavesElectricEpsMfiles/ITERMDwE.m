function [wE] = ITERMDwE(Einc,pG,E1data,E2data,input)
global nDIM it

   NS = input.NS;    
  AN1 = cell(1,NS);
  wEq = cell(2);      vEq = cell(2);      Eq = cell(2);
 g1Eq = cell(2);      rMq = cell(2);     rDq = cell(2);
 g1E  = cell(2,NS);    vE = cell(2,NS);      

itmax = 256;  
eta_M = 1 / (norm(E1data(:))^2 + norm(E2data(:))^2);  % normalization
   
% (1) Initial estimate of contrast sources
  [wE,E,r_M,ErrorM]=InitialEstimatewE(eta_M,Einc,pG,E1data,E2data,input);
  it = 0;   ErrorD = 1;
  disp(['Iteration = ', num2str(it)]);
  [CHI,r_D,eta_D,ErrorD] = UpdateEMContrast(wE,E,Einc,ErrorM,ErrorD,input);
  saveCHI(CHI);
  disp('--------------------------------------------------------------');
      
% (2) Iterative updating of contrast sources and permittivity contrast                 
  it = 1;  
  while (it <= itmax) 
    Norm_M = 0;
    for q = 1 : NS   
       AN1q = AN1{q};            
       for n = 1:nDIM 
         wEq{n} =  wE{n,q};   vEq{n} =  vE{n,q};     Eq{n} =   E{n,q};     
        g1Eq{n} = g1E{n,q};   rMq{n} = r_M{n,q};    rDq{n} = r_D{n,q};  
       end
  [wEq,vEq,Eq,g1Eq,AN1q,rMq] ...
      = InvITERwE(CHI,pG,wEq,vEq,Eq,g1Eq,AN1q,rMq,rDq,eta_M,eta_D,input);
      AN1{q} = AN1q;            
     for n = 1:nDIM 
          wE{n,q} =  wEq{n};    vE{n,q} = vEq{n};    E{n,q} =  Eq{n};   
         g1E{n,q} = g1Eq{n};   r_M{n,q} = rMq{n};        
     end
     Norm_M = Norm_M + norm(r_M{1,q}(:))^2 + norm(r_M{2,q}(:))^2;
   end %q_loop
     ErrorM = eta_M * Norm_M;
  
  disp(['Iteration = ', num2str(it)]);
  [CHI,r_D,eta_D,ErrorD]=UpdateEMContrast(wE,E,Einc,ErrorM,ErrorD,input);
  saveCHI(CHI);
  it = it+1; 
 end % while
 disp('--------------------------------------------------------------');
 Error  = sqrt(ErrorM + ErrorD);  
 disp(['Number of iterations = ' num2str(it)]);
 disp(['Total Error in M and D = ' num2str(Error)]);
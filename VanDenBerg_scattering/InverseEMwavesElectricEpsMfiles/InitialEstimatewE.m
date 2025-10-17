function [wE,E,r_M,ErrorM] = ...
                      InitialEstimatewE(eta_M,Einc,pG,E1data,E2data,input)
 NS = input.NS;   
 wE  = cell(2,NS);   E  = cell(2,NS);   r_M = cell(2,NS);  

% Determine contrast sources by back-projection and minimization of Norm_M
  Norm_M = 0; 
  for q = 1 : NS
    dummyq{1} = E1data(:,q);
    dummyq{2} = E2data(:,q); 
       g_M    = AdjDOPwE(pG,dummyq,input);
       A      = norm(g_M{1}(:))^2 + norm(g_M{2}(:))^2;
       GRv    = DOPwE(pG,g_M,input);  
       B      = norm(GRv{1}(:))^2 + norm(GRv{2}(:))^2; 
       alpha  = real(A/B);
     r_M{1,q} = E1data(:,q) - alpha * GRv{1};
     r_M{2,q} = E2data(:,q) - alpha * GRv{2};
      wE{1,q} = alpha * g_M{1}; 
      wE{2,q} = alpha * g_M{2};
    dummyq{1} = wE{1,q};
    dummyq{2} = wE{2,q};
       KwE    = KopE(dummyq,input);
      Norm_M  = Norm_M + norm(r_M{1,q}(:))^2 + norm(r_M{2,q}(:))^2;  
       E{1,q} = Einc{1,q} + KwE{1}; 
       E{2,q} = Einc{2,q} + KwE{2};
  end %q_loop
  
ErrorM = sqrt(eta_M * Norm_M);



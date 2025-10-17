function [dw] = ITERCGMDdw(data,input) 
global nDIM it 
NS    = input.NS;
AN1   = cell(1,NS);      
dwq   = cell(nDIM);      Pq = cell(nDIM);   vdwq = cell(nDIM);
g1dwq = cell(nDIM);     drq = cell(nDIM);  
g1dw  = cell(nDIM,NS);  vdw = cell(nDIM,NS); 

itmax = 256;   
eta_M = 1 / norm(data(:))^2;     % normalization factor eta_M

% (1) Initial estimate of interface sources and mass-density jumps
  [dw,P,P_inc,r_M,ErrorM] = InitialEstimatedw(eta_M,data,input);
   it = 0;  
   disp(['Iteration = ', num2str(it)]);                
  [Rfl,dr,eta_D,ErrorD] = UpdateContrastdw(dw,P,P_inc,ErrorM,input);
  saveRfl(Rfl);
  disp('--------------------------------------------------------------');
  
% (2) Iterative updates of interface sources and mass-density jumps
  it = 1;  
  while (it <= itmax)               
    Norm_M = 0;
    for q = 1 : NS   
      AN1q = AN1{q};            
      for n = 1:nDIM 
         dwq{n} = dw{n,q};       Pq{n} = P{n,q};    
       g1dwq{n} = g1dw{n,q};   vdwq{n} = vdw{n,q};     rMq = r_M{q};
         drq{n} = dr{n,q};
      end
      [dwq,vdwq,Pq,g1dwq,AN1q,rMq] ...
        = InvITERdw(Rfl,dwq,vdwq,Pq,g1dwq,AN1q,rMq,drq,eta_M,eta_D,input);
      AN1{q} = AN1q;
      for n=1:nDIM
         dw{n,q} = dwq{n};       P{n,q} = Pq{n};  
       g1dw{n,q} = g1dwq{n};   vdw{n,q} = vdwq{n};     r_M{q} = rMq;
      end  
      Norm_M = Norm_M + norm(r_M{q})^2; 
   end % qloop  
   ErrorM = eta_M * Norm_M;
   disp(['Iteration = ', num2str(it)]);       
   [Rfl,dr,eta_D,ErrorD] = UpdateContrastdw(dw,P,P_inc,ErrorM,input);   
   saveRfl(Rfl); 
   it = it+1;
 end % while
 
 disp('--------------------------------------------------------------');
 Error = sqrt(ErrorM + ErrorD); 
 disp(['Number of Iterations = ' num2str(it)]);
 disp(['Total Error in M and D = ' num2str(Error)]);
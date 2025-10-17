function [w,CHI] = ITERCGMDreg(u_inc,G_R,data,input)
global it
FFTG = input.FFTG;  NS = input.NS;    
v = cell(1,NS);    g_1 = cell(1,NS);    AN_1 = cell(1,NS);

itmax = 256;  
eta_M = 1 / norm(data(:))^2;     % normalization factor eta_M

[w,u,CHI,r_M,r_D,eta_D,ErrorD,ErrorM] = ...
                           InitialEstimate(eta_M,u_inc,G_R,data,input);            
it = 1;  
while (it <= itmax)         % alternate updating of w and CHI   
 Norm_M = 0;
 for q = 1 : NS
   g  = eta_M * AdjDopM(G_R,r_M{q},input)  ...
         + eta_D *(r_D{q}-AdjKop(conj(CHI).*r_D{q},FFTG));  
      g = (abs(CHI)>=0.001) .*  g; 
   AN = norm(g(:))^2;
   if it > 1 ; BN = sum( g(:) .* conj(g_1{q}(:)) ); end 
   if it == 1  
      v{q} = g; 
   else     
      v{q} = g + real((AN-BN)/AN_1{q}) * v{q};  
   end 
     g_1{q} = g; 
    AN_1{q} = AN;   
  % determine step length alpha    
  AN    = sum( g(:) .* conj(v{q}(:)) ); 
  GRv   = DopM(G_R,v{q},input);
  BN_M  = norm(GRv(:))^2;  
  Kv    = Kop(v{q},FFTG);
  Lv    = v{q} - CHI.* Kv;
  BN_D  = norm(Lv(:))^2;
  alpha = real( AN / (eta_M * BN_M + eta_D * BN_D) );   
  % update contrast source w and  residual  error  in M
    w{q} =   w{q} + alpha *   v{q};
    u{q} =   u{q} + alpha *  Kv; 
  r_M{q} = r_M{q} - alpha * GRv(:);  
  Norm_M = Norm_M + norm(r_M{q}(:))^2; 
 end %q_loop
 ErrorM = eta_M * Norm_M;
  
 % Update contrast  by  minimization of  Norm_D = || CHI * u - w ||^2 
  disp(['Iteration = ', num2str(it)]);                   saveCHI(CHI);
  [CHI,r_D,eta_D,ErrorD]=UpdateContrastReg(w,u,u_inc,ErrorM,ErrorD,input);
  disp('--------------------------------------------------------------');
  it = it+1; 
end % while

Error  = sqrt(ErrorM + ErrorD);  
disp(['Number of iterations = ' num2str(it)]);
disp(['Total Error in M and D = ' num2str(Error)]);
function [wE,vE,E,gE,AN,r_M] ...
       = InvITERwE(CHI,pG,wE,vE,E,g_1E,AN_1,r_M,r_D,eta_M,eta_D,input)       
global it  
% determine conjugate gradient directions 
    g_M    = AdjDOPwE(pG,r_M,input);
  dummy{1} = conj(CHI) .* r_D{1};
  dummy{2} = conj(CHI) .* r_D{2}; 
   [adjKE] = AdjKopE(dummy,input);
    gE{1}  = eta_M * g_M{1} + eta_D * (r_D{1} - adjKE{1}); 
    gE{2}  = eta_M * g_M{2} + eta_D * (r_D{2} - adjKE{2}); 
    AN     = norm(gE{1}(:))^2 + norm(gE{2}(:))^2;
if it > 1 
    BN  = sum(gE{1}(:).*conj(g_1E{1}(:))+gE{2}(:).*conj(g_1E{2}(:)));
end 
if it == 1  
    vE = gE; 
else     
    vE{1} = gE{1} + real((AN-BN)/AN_1) * vE{1};  
    vE{2} = gE{2} + real((AN-BN)/AN_1) * vE{2};  
end 
% determine step length alpha    
  AN_D  = sum(gE{1}(:).*conj(vE{1}(:))+gE{2}(:).*conj(vE{2}(:)));
  Gv_M  = DOPwE(pG,vE,input);
  BN_M  = norm(Gv_M{1}(:))^2 + norm(Gv_M{2}(:))^2;  
  [KvE] = KopE(vE,input);
 LvE{1} = vE{1} - CHI .* KvE{1};
 LvE{2} = vE{2} - CHI .* KvE{2};
  BN_D  = norm(LvE{1}(:))^2 + norm(LvE{2}(:))^2;
  alpha = real(AN_D) / (eta_M * BN_M + eta_D * BN_D);  
% update contrast source wE and field E
  wE{1} = wE{1}  + alpha * vE{1};
  wE{2} = wE{2}  + alpha * vE{2};
   E{1} = E{1}   + alpha * KvE{1}; 
   E{2} = E{2}   + alpha * KvE{2};
% update residual error in M    
 r_M{1} = r_M{1} - alpha * Gv_M{1}(:);
 r_M{2} = r_M{2} - alpha * Gv_M{2}(:);
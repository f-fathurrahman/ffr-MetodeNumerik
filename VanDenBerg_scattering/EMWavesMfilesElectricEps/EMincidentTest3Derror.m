clear all; clc; close all; clear workspace
input = initEM();  gam0 = input.gamma_0;  xS = input.xS;
    a = input.a;  theta = 0:.01:pi;      phi = 0:.01:2*pi; 
    [THETA,PHI] = ndgrid(theta,phi);
% (1) Cartesian coordinates: field from electric dipole in negative X_1 --
   X1   = a*sin(THETA).*cos(PHI)-xS(1); X2 = a*sin(THETA).*sin(PHI)-xS(2);    
   X3   = a*cos(THETA)-xS(3); 
  DIS   = sqrt(X1.^2 + X2.^2 + X3.^2);
   X1   = X1 ./ DIS;        X2 = X2 ./ DIS;         X3 = X3 ./ DIS;
    G   = exp(-gam0*DIS) ./ (4*pi*DIS);
    dG  = -(gam0 + 1./DIS) .* G;
   dG11 = ((3*X1.*X1 -1).* (gam0./DIS + 1./DIS.^2) + gam0^2 *X1.*X1) .* G;
   dG21 = ( 3*X2.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2 *X2.*X1) .* G;
   dG31 = ( 3*X3.*X1    .* (gam0./DIS + 1./DIS.^2) + gam0^2 *X3.*X1) .* G;                            
  E{1}  = - (-gam0^2 * G + dG11);
  E{2}  = - dG21;
  E{3}  = - dG31;
  ZH{2} = -gam0 * X3 .* dG;
  ZH{3} =  gam0 * X2 .* dG;

[E1,E2,E3,ZH2,ZH3] = IncEMsphereError(input,THETA,PHI);
  
  % Print the normalized error on the sphere boundary at r = a
  norm_Error = norm(E{1}(:)-E1(:),1) / norm(E{1}(:),1);      
               disp(['Error ',' E_1',' = ',num2str(norm_Error)]); 
  norm_Error = norm(E{2}(:)-E2(:),1) / norm(E{2}(:),1);      
               disp(['Error ',' E_2',' = ',num2str(norm_Error)]); 
  norm_Error = norm(E{3}(:)-E3(:),1) / norm(E{3}(:),1);      
               disp(['Error ',' E_3',' = ',num2str(norm_Error)]);              
                  
  norm_Error = norm(ZH{2}(:)-ZH2(:),1) / norm(ZH{2}(:),1);      
               disp(['Error ','ZH_2',' = ',num2str(norm_Error)]);             
  norm_Error = norm(ZH{3}(:)-ZH3(:),1) / norm(ZH{3}(:),1);      
               disp(['Error ','ZH_3',' = ',num2str(norm_Error)]); 
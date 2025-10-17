function [input] = AcousticModifiedContrast(input)
global nDIM; 

if     nDIM == 1;   R = sqrt(input.X1.^2);   
elseif nDIM == 2;   R = sqrt(input.X1.^2 + input.X2.^2);
elseif nDIM == 3;   R = sqrt(input.X1.^2 + input.X2.^2 + input.X3.^2);
end % if

% (1) Compute wave speed contrast CHI^c -----------------------------------  
   ratio_c_sct2 = (1 + (input.c_sct.^2/input.c_0^2 - 1) * (R < input.a));
   input.CHI_c  = 1 - (1./ratio_c_sct2);
   
% (2) Compute gradient mass density contrast CHI^\rho/drho----------------- 
   ratio_rho_sct =  1 + (input.rho_sct/input.rho_0-1) * (R < input.a);
   inverse_ratio_rho   =  1./ ratio_rho_sct;
              drho_sct = GRAD(inverse_ratio_rho, input);
   for n =1:nDIM
     input.CHI_drho{n} = ratio_rho_sct.* drho_sct{n};
   end
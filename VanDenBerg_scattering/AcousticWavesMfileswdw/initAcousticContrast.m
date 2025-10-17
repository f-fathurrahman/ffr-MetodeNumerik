function [input] = initAcousticContrast(input)
global nDIM; 

input.a  = 40;  % half width slab / radius circle cylinder / radius sphere

if     nDIM == 1;   R = sqrt(input.X1.^2);   
elseif nDIM == 2;   R = sqrt(input.X1.^2 + input.X2.^2);
elseif nDIM == 3;   R = sqrt(input.X1.^2 + input.X2.^2 + input.X3.^2);
end % if

% (1) Compute compressibbily contrast (kappa = 1/(rho*c^2) ---------------
kappa_0        = 1 / (input.rho_0   * input.c_0^2);
kappa_sct      = 1 / (input.rho_sct * input.c_sct^2);
contrast_kappa = 1 - kappa_sct /kappa_0;
input.CHI_kap  = contrast_kappa * (R < input.a);   

% (2) Compute mass density contrast --------------------------------------
contrast_rho   = 1 - input.rho_sct / input.rho_0; 
input.CHI_rho  = contrast_rho * (R < input.a);
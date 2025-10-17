function [input] = initEMcontrast(input)
global nDIM; 
input.a  = 40;      % radius circle cylinder / radius sphere 

if     nDIM == 2;   R = sqrt(input.X1.^2 + input.X2.^2);
elseif nDIM == 3;   R = sqrt(input.X1.^2 + input.X2.^2 + input.X3.^2);
end % if

% (1) Compute permittivity contrast --------------------------------------
input.CHI_eps = (1-input.eps_sct) * (R < input.a);  

% (2) Compute permeability contrast --------------------------------------
input.CHI_mu  = (1-input.mu_sct)  * (R < input.a);  
function [Ercv] = DOPwE(pG,wE,input)

gam0 = input.gamma_0;    dx = input.dx;  

Ercv{1} = zeros(input.NR,1); 
Ercv{2} = zeros(input.NR,1);
for p = 1 : input.NR     
 Ercv{1}(p,1) = dx^2*sum((gam0^2*pG.GR{p}(:)-pG.dGR11{p}(:)).* wE{1}(:) ...
                                           -pG.dGR21{p}(:) .* wE{2}(:)); 
 Ercv{2}(p,1) = dx^2*sum((gam0^2*pG.GR{p}(:)-pG.dGR22{p}(:)).* wE{2}(:) ...
                                           -pG.dGR21{p}(:) .* wE{1}(:));
end; % p_loop

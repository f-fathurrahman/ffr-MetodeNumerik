function [GMw] = DOPMwZv(dG_R,w_Zv,input)

gam0 = input.gamma_0;    
dx   = input.dx;   

GMw = zeros(input.NR,1);
for p = 1 : input.NR    
    GMw(p,1) =  gam0 * dx^2 * sum(dG_R{1,p}(:) .* w_Zv{1}(:)) ...
              + gam0 * dx^2 * sum(dG_R{2,p}(:) .* w_Zv{2}(:));       
end % p_loop 
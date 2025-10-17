function [adjGMZv] = AdjDOPMwZv(dG_R,f,input)
 
gam0 = input.gamma_0;    
dx   = input.dx;

adjGMZv = cell(1,2);  
adjGMZv{1} = zeros(input.N1,input.N2);  
adjGMZv{2} = zeros(input.N1,input.N2);
for p = 1: input.NR  
    adjGMZv{1} = adjGMZv{1} + conj(gam0 * dx^2 * dG_R{1,p}) * f(p); 
    adjGMZv{2} = adjGMZv{2} + conj(gam0 * dx^2 * dG_R{2,p}) * f(p);
end % p_loop                                                
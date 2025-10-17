function [adjGR] = AdjDopM(G_R,f,input)

gam0 = input.gamma_0;    
dx   = input.dx;
 
adjGR = zeros(input.N1,input.N2);     
for p = 1: input.NR
    adjGR = adjGR + conj(gam0^2 * dx^2 *  G_R{p}) * f(p);     
end % p_loop
                                                
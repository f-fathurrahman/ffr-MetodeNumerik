function [GRw] = DopM(G_R,w,input)

gam0 = input.gamma_0;    
dx   = input.dx;
 
GRw = zeros(input.NR,1);
for p = 1: input.NR; 
     GRw(p,1) = (gam0^2 * dx^2) * sum(G_R{p}(:).*w(:));     
end; % p_loop 
                                             
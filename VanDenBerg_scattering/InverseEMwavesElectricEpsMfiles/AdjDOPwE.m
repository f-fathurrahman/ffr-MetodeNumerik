function [adjGR] = AdjDOPwE(pG,f,input)

gam0 = input.gamma_0;    dx = input.dx;  
 
adjGR{1} = zeros(input.N1,input.N2);  
adjGR{2} = zeros(input.N1,input.N2);
for p = 1 : input.NR     
  adjGR{1} = adjGR{1} + conj(gam0^2*pG.GR{p}-pG.dGR11{p}) .* f{1}(p) ...
                                      - conj(pG.dGR21{p}) .* f{2}(p); 
  adjGR{2} = adjGR{2} + conj(gam0^2*pG.GR{p}-pG.dGR22{p}) .* f{2}(p) ...
                                      - conj(pG.dGR21{p}) .* f{1}(p);
end % p_loop
adjGR{1} = dx^2 * adjGR{1};
adjGR{2} = dx^2 * adjGR{2};
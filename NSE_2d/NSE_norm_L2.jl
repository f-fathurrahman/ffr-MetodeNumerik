# Computes the L2 norm on a 2D domain
# ||u||_2 = (\int_Omega  u^2 dx dy)^{1/2}
function NSE_norm_L2(u, dx, dy)
    #global dx dy;
    return sqrt(sum(u .* u) *dx*dy)
end
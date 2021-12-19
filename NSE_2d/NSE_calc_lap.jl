# Computes the discrete Laplacian ∇²u
function NSE_calc_lap(u, dx, dy, imm, ip, jp, jm, ic, jc)
    #global dx dy
    #global im ip jp jm ic jc
    hc = ( u[ip,jc] - 2*u[:,:] + u[imm,jc] )/(dx*dx) + ( u[ic,jp] - 2*u[:,:] + u[ic,jm] )/(dy*dy)
    return hc
end

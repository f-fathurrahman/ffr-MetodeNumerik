"""
meshgrid(vx,vy)
Computes an (x,y)-grid from the vectors (vx,vy).
For more information, see the MATLAB documentation.

From: https://github.com/ChrisRackauckas/VectorizedRoutines.jl
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where {T}
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    return (repeat(vx, m, 1), repeat(vy, 1, n))
end


# Initial condition for a 2D periodic jet
# study of the Kelvin-helmholtz instability
function NSE_F_init_KH!(
    region,
    Lx, Ly, x, y,
    U0, Rj, Pj, Ax, lamx,
    u
)
    
    dx = region.dx
    dy = region.dy

    idx_im = region.idx_im
    idx_ip = region.idx_ip
    idx_jp = region.idx_jp
    idx_jm = region.idx_jm
    idx_ic = region.idx_ic
    idx_jc = region.idx_jc
    
    # 2D grid
    # XXX FF: why compute this again?
    #     Probably because the grid can be mixed, e.g xc with ym
    xx, yy = meshgrid(x, y);
    xx = xx'
    yy = yy'
    
    # local y coordinate
    rr = yy .- Ly/2;
    
    Nxm = region.Nxm
    Nym = region.Nym
    # velocity u
    for j in 1:Nym, i in 1:Nxm
        u[i,j] = U0 * 0.5 * ( 1.0 + tanh(0.5*Pj*(1 - abs(rr[i,j])/Rj)))
        #perturbation in $x$ direction
        u[i,j] += Ax * u[i,j] * sin(2*pi/lamx * xx[i,j])
    end

    return
end


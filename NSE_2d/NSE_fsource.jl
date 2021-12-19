function NSE_fsource(Lx, Ly, x , y)
     aa = 2π/Lx
     bb = 2π/Ly
     fs = (aa*aa + bb*bb) * ( cos.(aa*x) .* sin.(bb*y) )
     return fs
end

# -*- coding: utf-8 -*-
import Plots, PlotThemes
Plots.theme(:dark)

function propagate_Ex_Hy(Ntime::Int64)
    Nz = 200
    # Initially, there is no electric nor magnetic field
    Ex = zeros(Float64, Nz)
    Hy = zeros(Float64, Nz)
    # Pulse parameters
    idx_middle_z = Int64(Nz/2)
    t0 = 40.0
    spread = 12.0;
    # Main FDTD Loop
    for n in 1:Ntime
        # Calculate the Ex field
        for k in 2:Nz
            Ex[k] = Ex[k] + 0.5 * (Hy[k-1] - Hy[k])
        end
        # Put a Gaussian pulse in the middle
        # This pulse only depends on time, no spatial dependence (it is fixed at approximately x=0)
        pulse = exp(-0.5*( (t0 - n) / spread)^2) # Δt is implicitly 1
        Ex[idx_middle_z] = pulse
        # Calculate the Hy field
        for k in 1:(Nz-1)
            Hy[k] = Hy[k] + 0.5 * (Ex[k] - Ex[k + 1])
        end
    end
    return Ex, Hy
end;

Ex, Hy = propagate_Ex_Hy(10)
f = Plots.plot(Ex, fmt=:svg, size=(500,300), label="Ex")
Plots.plot!(Hy, label="Hy")
Plots.ylims!(f, (-1.0, 1.0))


Ex, Hy = propagate_Ex_Hy(50)
f = Plots.plot(Ex, fmt=:svg, size=(500,300), label="Ex")
Plots.plot!(Hy, label="Hy")
Plots.ylims!(f, (-1.0, 1.0))

Ex, Hy = propagate_Ex_Hy(75)
f = Plots.plot(Ex, fmt=:svg, size=(500,300), label="Ex")
Plots.plot!(Hy, label="Hy")
Plots.ylims!(f, (-1.0, 1.0))

Ex, Hy = propagate_Ex_Hy(300)
f = Plots.plot(Ex, fmt=:svg, size=(500,300), label="Ex")
Plots.plot!(Hy, label="Hy")
Plots.ylims!(f, (-1.0, 1.0))

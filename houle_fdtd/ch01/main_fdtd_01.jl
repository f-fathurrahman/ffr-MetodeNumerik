
import Plots, PlotThemes
Plots.theme(:dark)

# This is necessary for GR no-display mode (for faster writing to file)
ENV["GKSwstype"] = "nul"

# evaluate a time dependent pulse, position will be decided outside this function
function evaluate_pulse(t; spread=12.0, t0=40.0)
    return exp(-0.5*( (t0 - t) / spread)^2)
end


function main_fdtd()
    Nz = 200
    # Initially, there is no electric nor magnetic field
    Ex = zeros(Float64, Nz)
    Hy = zeros(Float64, Nz)
    # Pulse parameters
    idx_middle_z = Int64(Nz/2)
    t0 = 40.0
    spread = 12.0
    Δt = 1.0
    Ntime = 100
    # Main FDTD Loop
    for n in 1:Ntime
        # Calculate the Ex field
        for k in 2:Nz
            Ex[k] = Ex[k] + 0.5 * (Hy[k-1] - Hy[k])
        end
        pulse = evaluate_pulse(n*Δt)
        Ex[idx_middle_z] += pulse # ffr: add the pulse instead of modify it
        #
        # Calculate the Hy field
        for k in 1:(Nz-1)
            Hy[k] = Hy[k] + 0.5 * (Ex[k] - Ex[k + 1])
        end
        fig = Plots.plot(Ex, size=(1000,500))
        Plots.plot!([idx_middle_z], [pulse], marker=:o)
        Plots.plot!(fig, Hy)
        Plots.ylims!(fig, (-1.0, 1.0))
        Plots.savefig(fig, "IMG_$(n).png")
        println("Done $n")
    end
end


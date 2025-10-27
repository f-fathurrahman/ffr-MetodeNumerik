# -*- coding: utf-8 -*-
import Plots, PlotThemes
Plots.theme(:dark)

p = Plots.plot([sin, cos], zeros(0), leg = false, xlims = (0, 2π), ylims = (-1, 1))
anim = Plots.Animation()
for x = range(0, stop = 2π, length = 20)
    # add new figure
    push!(p, x, Float64[sin(x), cos(x)])
    Plots.frame(anim)
end

Plots.gif(anim, "anim_fps15.gif", fps = 15)

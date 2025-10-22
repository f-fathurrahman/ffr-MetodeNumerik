# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.5.1
#     language: julia
#     name: julia-1.5
# ---

# # Second Order
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Vec

p = let
	d = 0.2
	e = 0.02
	f = x -> x^2 + d*x^3 - e*x^4
	f′ = x -> 2x + 3d*x^2 - 4e*x^3
	f″ = x -> 2 + 6d*x - 12e*x^2
	xdomain = (-2,2)

	x = 1.25

	G = GroupPlot(2, 1, groupStyle="horizontal sep=1.25cm", style="ymin=-1, xmin=-2.25, width=6cm, ymax=$(f(2)), axis lines=left, xtick=\\empty, ytick=\\empty, xlabel=\$x\$, ylabel=\$f\$") # $

	g = α -> f(x) + (α-x)*f′(x)
	plots = Plots.Plot[]
	push!(plots, Plots.Linear(f, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x], [f(x)], style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	push!(G, Axis(plots))

	g = α -> f(x) + (α-x)*f′(x) + (α-x)^2*f″(x)/2
	plots = Plots.Plot[]
	push!(plots, Plots.Linear(f, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x], [f(x)], style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	push!(G, Axis(plots))

	G
end

plot(p)

# +
using Vec

p = let
	d = 0.2
	e = 0.02
	f = x -> x^2 + d*x^3 - e*x^4
	f′ = x -> 2x + 3d*x^2 - 4e*x^3
	f″ = x -> 2 + 6d*x - 12e*x^2
	xdomain = (-2,2)

	x = 1.25
	x2 = x - f′(x)/f″(x)

	G = GroupPlot(2, 1, groupStyle="horizontal sep=1.25cm", style="ymin=-1, xmin=-2.25, width=6cm, ymax=$(f(2)), axis lines=left, xtick={$x, $x2}, xticklabels={\$x^{(k)}\$, \$x^{(k+1)}\$}, xlabel=\$x\$, xticklabel style={text height=2ex}") # $

	g = α -> f(x) + (α-x)*f′(x) + (α-x)^2*f″(x)/2
	plots = Plots.Plot[]
	push!(plots, Plots.Linear(f, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x, x2], [f(x), f(x2)], style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	push!(G, Axis(plots, ylabel=L"f", style="ytick=\\empty"))

	g = α -> f′(x) + (α-x)*f″(x)
	plots = Plots.Plot[]
	push!(plots, Plots.Linear([xdomain[1], xdomain[2]], [0,0], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear([x2,x2], [-2,6], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear(f′, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x, x2], [f′(x), f′(x2)], style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	push!(G, Axis(plots, ylabel=L"f'", xmin=-2, ymin=-2, ymax=6, style="ytick={0}"))

	G
end

plot(p)

# +
p = let

	G = GroupPlot(3, 1, groupStyle="horizontal sep=1.25cm, xlabels at=edge bottom, ylabels at=edge left", style="xlabel=\$x\$, ylabel=\$f\$, width=6cm, axis lines=left, xticklabel style={text height=2ex}")

	f = x -> log(exp(x)+1) - x/2 + 1
	f′ = x -> 1/(1+exp(-x)) - 0.5
	f″ = x -> exp(x)/(exp(x)+1)^2

	x = 2.17732
	x2 = x - f′(x)/f″(x)
	xdomain = (-4,4)
	ydomain = (1.15,3.25)

	g = α -> f(x) + (α-x)*f′(x) + (α-x)^2*f″(x)/2
	g2 = α -> f(x2) + (α-x2)*f′(x2) + (α-x2)^2*f″(x2)/2
	plots = Plots.Plot[]
	push!(plots, Plots.Linear([x,x], [ydomain[1], ydomain[2]], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear([x2,x2], [ydomain[1], ydomain[2]], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear(f, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x,x2], f.([x,x2]), style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	push!(plots, Plots.Linear(g2, xdomain, style="solid, thick, pastelPurple, mark=none"))
	push!(G, Axis(plots, title="Oscillation", xmin=xdomain[1], xmax=xdomain[2], ymin=ydomain[1], ymax=ydomain[2], style="xtick=\\empty, ytick=\\empty"))

	f = x -> x < -π/2 ? sin(x) : 2sin((x-π/2)/2) + 1
	f′ = x -> cos(x)
	f″ = x -> -sin(x)

	x = -π + 0.2
	x2 = x - f′(x)/f″(x)
	xdomain = (-4.25,3)
	ydomain = (-3,2)

	g = α -> f(x) + (α-x)*f′(x) + (α-x)^2*f″(x)/2
	plots = Plots.Plot[]
	push!(plots, Plots.Linear([x,x], [ydomain[1], ydomain[2]], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear([x2,x2], [ydomain[1], ydomain[2]], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear(f, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x,x2], f.([x,x2]), style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	ax = Axis(plots, title="Overshoot", xmin=xdomain[1], xmax=xdomain[2], ymin=ydomain[1], ymax=ydomain[2],
	          style="xtick={$x, $(x2)}, xticklabels={\$x^{(k)}\$, \$x^{(k+1)}\$}, ytick=\\empty, ")
	push!(G, ax)

	c = -π - 0.4
	f = x -> x > c ? sin(x) : cos(c)*(x-c) + sin(c) - 0.1*abs(x-c)^2
	f′ = x -> cos(x)
	f″ = x -> -sin(x)

	x = -π - 0.3
	x2 = x - f′(x)/f″(x)
	xdomain = (-8,-1)
	ydomain = (-2,4)

	g = α -> f(x) + (α-x)*f′(x) + (α-x)^2*f″(x)/2
	plots = Plots.Plot[]
	push!(plots, Plots.Linear([x,x], [ydomain[1], ydomain[2]], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear([x2,x2], [ydomain[1], ydomain[2]], style="solid, gray, mark=none"))
	push!(plots, Plots.Linear(f, xdomain, style="solid, thick, black, mark=none"))
	push!(plots, Plots.Scatter([x,x2], f.([x,x2]), style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear(g, xdomain, style="solid, thick, pastelBlue, mark=none"))
	ax = Axis(plots, title="Negative \$f''\$", xmin=xdomain[1], xmax=xdomain[2], ymin=ydomain[1], ymax=ydomain[2],
	          style="xtick={$x, $(x2)}, xticklabels={\$x^{(k)}\$, \$x^{(k+1)}\$}, ytick=\\empty, ")
	push!(G, ax)

	G
end

plot(p)

# +
using Vec
using LinearAlgebra

function _line_search(f, x, d)
    d = normalize(d)
    objective = α -> f(x + α*d)
    v, α = f(x), 1e-6
    while f(x + α*d) < v
        v = f(x + α*d)
        α += 1e-6
    end
    return x + α*d
end

abstract type DescentMethod end
mutable struct DFP <: DescentMethod
	Q
end
function init!(M::DFP, f, ∇f, x)
	M.Q = Matrix(1.0I, length(x), length(x))
	return M
end
function step!(M::DFP, f, ∇f, x)
	Q, g = M.Q, ∇f(x)
	x′ = _line_search(f, x, -Q*g)
	g′ = ∇f(x′)
	δ = x′ - x
    γ = g′ - g
    Q[:] = Q - Q*γ*γ'*Q/(γ'*Q*γ) + δ*δ'/(δ'*γ)
    return x′
end

mutable struct myBFGS <: DescentMethod
	Q
end
function init!(M::myBFGS, f, ∇f, x)
	M.Q = Matrix(1.0I, length(x), length(x))
	return M
end
function step!(M::myBFGS, f, ∇f, x)
	Q, g = M.Q, ∇f(x)
	x′ = _line_search(f, x, -Q*g)
	g′ = ∇f(x′)
	δ = x′ - x
    γ = g′ - g
    Q[:] = Q - (δ*γ'*Q + Q*γ*δ')/(δ'*γ) + (1 + (γ'*Q*γ)/(δ'*γ))[1]*(δ*δ')/(δ'*γ)
    return x′
end

mutable struct LimitedMemoryBFGS <: DescentMethod
	m
	δs
	γs
	qs
end
function init!(M::LimitedMemoryBFGS, f, ∇f, x)
	M.δs = []
	M.γs = []
    M.qs = []
	return M
end
function step!(M::LimitedMemoryBFGS, f, ∇f, x)
    δs, γs, qs, g = M.δs, M.γs, M.qs, ∇f(x)
    m = length(δs)
    if m > 0
        q = g
        for i in m : -1 : 1
            qs[i] = copy(q)
            q -= (δs[i]⋅q)/(γs[i]⋅δs[i])*γs[i]
        end
        z = (γs[m] .* δs[m] .* q) / (γs[m]⋅γs[m])
        for i in 1 : m
            z += δs[i]*(δs[i]⋅qs[i] - γs[i]⋅z)/(γs[i]⋅δs[i])
        end
        x′ = _line_search(f, x, -z)
    else
        x′ = _line_search(f, x, -g)
    end
    g′ = ∇f(x′)
    push!(δs, x′ - x); push!(γs, g′ - g)
    push!(qs, zeros(length(x)))
    while length(δs) > M.m
        popfirst!(δs); popfirst!(γs); popfirst!(qs)
    end
    return x′
end

p = let
	f = x -> (1-x[1])^2 + 5*(4x[2] - x[1]^2)^2
	∇f = x -> [2*(10x[1]^3 - 40x[1]*x[2] + x[1] - 1), -40*(x[1]^2 - 4x[2])]

	xdomain = (-3, 2)
	ydomain = (-0.5, 2)

	function this_step!(M::DescentMethod, v::VecE2{Float64})
	    x = convert(Vector{Float64}, v)
	    return VecE2{Float64}(step!(M, f, ∇f, x))
	end
	function run_descent_method(M::DescentMethod, x₀::VecE2{Float64}, N::Int)
	    pts = [x₀]
	    init!(M, f, ∇f, convert(Vector{Float64}, x₀))
	    for i in 1 : N
	        push!(pts, this_step!(M, pts[end]))
	    end
	    return pts
	end
	function get_descent_plot(pts::Vector{VecE2{Float64}}, name, color::String="black")
	    Plots.Linear([P.x for P in pts], [P.y for P in pts], style="thick, $color, solid, mark=none, line join=round", legendentry=name) #$
	end

	x₀ = VecE2{Float64}(-1,1.75)
	N = 15
	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[2,10,50,200,500], style="forget plot, width=\\columnwidth", xbins=150, ybins=150))

	stuff = Tuple{DescentMethod, String, String}[]
	push!(stuff, (DFP(NaN), "DFP", "black"))
	push!(stuff, (myBFGS(NaN), "BFGS", "gray"))
	push!(stuff, (LimitedMemoryBFGS(3,NaN,NaN,NaN), "L-BFGS (\$m=3\$)", "pastelGreen"))
	push!(stuff, (LimitedMemoryBFGS(2,NaN,NaN,NaN), "L-BFGS (\$m=2\$)", "pastelSeaGreen"))
	push!(stuff, (LimitedMemoryBFGS(1,NaN,NaN,NaN), "L-BFGS (\$m=1\$)", "pastelBlue"))
	for (M, name, color) in stuff
	    pts = run_descent_method(M, x₀, N)
	    push!(plots, get_descent_plot(pts, name, color))
	end
	push!(plots, Plots.Scatter([1], [1/4], style="mark=*, mark size = 1, mark options={draw=black, fill=black}"))
	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, view={0}{90}, legend pos=outer north east")
end
plot(p)

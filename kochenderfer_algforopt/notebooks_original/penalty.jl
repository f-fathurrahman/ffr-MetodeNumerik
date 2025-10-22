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

# # Penalty
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
p = let

	f = x -> x[1]^2 - 0.9*x[1]*x[2] + x[2]^2

	xdomain = (-2, 1)
	ydomain = (-2, 1)


	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[0.1,0.2,0.5,1,2,3,4], style="width=\\columnwidth"))
	push!(plots, Plots.Linear3([-1.5,-1.5], [-2.0, 1.0], [5,5], style="black, dotted, mark=none"))
	push!(plots, Plots.Linear3([-1.5,-1.5], [-1.5, 0.5], [6,6], style="black, solid, mark=none"))
	push!(plots, Plots.Linear3([ 0.5, 0.5], [-2.0, 1.0], [5,5], style="black, dotted, mark=none"))
	push!(plots, Plots.Linear3([ 0.5, 0.5], [-1.5, 0.5], [6,6], style="black, solid, mark=none"))
	push!(plots, Plots.Linear3([-2.0, 1.0], [-1.5,-1.5], [5,5], style="black, dotted, mark=none"))
	push!(plots, Plots.Linear3([-1.5, 0.5], [-1.5,-1.5], [6,6], style="black, solid, mark=none"))
	push!(plots, Plots.Linear3([-2.0, 1.0], [ 0.5, 0.5], [5,5], style="black, dotted, mark=none"))
	push!(plots, Plots.Linear3([-1.5, 0.5], [ 0.5, 0.5], [6,6], style="black, solid, mark=none"))

	Axis(plots, width="1.1*9cm", height="1.1*9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +
p = let

	f = x -> begin
		r = norm(x)
		θ = atan(x[2], x[1])
		(sin(3θ) + r - 1)
	end

	xdomain = (-2,2)
	ydomain = (-2,2)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[-1, -0.5, 0.5, 1, 1.5, 2], labels=false))
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[0], contour_style="draw color=black"))

	Axis(plots, width="1.1*9cm", height="1.1*9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, axis equal, view={0}{90}")
end

plot(p)

# +
using Optim

p = let

	f = x ->  -exp(-(x[1]*x[2] - 1.5)^2 -(x[2]-1.5)^2)
	h = x -> x[1] - x[2]^2

	xdomain = (0, 3)
	ydomain = (0, 3)

	x = optimize(x->f(x) + 10000*abs(h(x))^2, [1.0,1.0]).minimizer

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, xbins=101, ybins=101, levels=[-0.1,-0.2,-0.3,-0.4,-0.5,-0.7,f(x),-0.97]))
	push!(plots, Plots.Contour(h, xdomain, ydomain, xbins=101, ybins=101, levels=[0], contour_style="draw color=black"))
	push!(plots, Plots.Scatter([x[1]], [x[2]], style="black, mark size=1, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Node(L"h(\vect x) = 0", 2.75, sqrt(2.75), axis="axis cs", style="anchor=south east"))
	push!(plots, Plots.Node(L"\vect x^*",x[1],x[2], axis="axis cs", style="anchor=north west"))

	Axis(plots, xmin=xdomain[1], xmax=xdomain[2], ymin=ydomain[1], ymax=ydomain[2],
	     width="9cm", height="9cm", xlabel=L"x_1", ylabel=L"x_2",
	     style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +
using Optim
import Printf: @sprintf
p = let
	f = x -> sin(x)
	# x^2 <= 1
	L = (x,μ) -> f(x) + μ*(x^2 - 1)
	dom = (-4,4)

	x_star = -1

	plots = Plots.Plot[]
	for (i,μ) in enumerate(collect(0:8)*0.25)
	    g = x->L(x,μ)
	    color = "pastelPurple!$(109-9i)"
	    push!(plots, Plots.Linear(g, dom, style="solid, thick, $color, mark=none", legendentry=@sprintf("\$\\mu = %.2f\$", μ)))
	    x = optimize(g, -4, 4).minimizer
	    push!(plots, Plots.Scatter([x], [g(x)], style="mark=*, mark size=1, mark options={draw=$color, fill=$color}, forget plot"))
	end
	push!(plots, Plots.Linear(f, dom, style="solid, black, mark=none"))
	push!(plots, Plots.Linear(f, (-1,1), style="solid, thick, pastelBlue, mark=none"))
	push!(plots, Plots.Scatter([x_star], [f(x_star)], style="mark=*, mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}"))
	push!(plots, Plots.Node(L"x^*", x_star, f(x_star), style="anchor=south west"))
	Axis(plots, style="xlabel=\$x\$, ylabel=\$y\$, enlarge x limits=0, ymax=10, legend cell align=left, legend style={draw=black, fill=white}, legend pos=outer north east, axis on top", width="7cm")
end

plot(p)

# +
using Optim

p = let

	G = GroupPlot(1,2,groupStyle="vertical sep=1cm")

	# PLOT 1

	f = x -> x[1] + x[2] + x[1]*x[2]
	h = x -> x[1]^2 + x[2]^2 - 1

	xdomain = (-2, 2)
	ydomain = (-2, 2)

	s = (√2-1)/(√2-2)
	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[f([-1,0]),0,f([(1+√2)/(2+√2),(1+√2)/(2+√2)]),f([s,s]), -2, 1, 3]))
	push!(plots, Plots.Contour(h, xdomain, ydomain, levels=[0], contour_style="draw color=black, labels=false"))
	push!(plots, Plots.Scatter([-1,0,(1+√2)/(2+√2),s], [0,-1,(1+√2)/(2+√2),s], style="black, mark size=1.5, mark options={draw=black, fill=black}"))
	push!(G, Axis(plots, width="9cm", height="9cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}"))

	# PLOT 2

	λ_opt = 1/2
	push!(G,
		Axis(
			[Plots.Linear(λ -> -λ - 1/(2*λ+1), (1/2,5), style="solid, black, mark=none"),
			 Plots.Linear([λ_opt,λ_opt],[-20,5], style="solid, black!40, mark=none"),
 				 Plots.Node(L"\lambda^*", λ_opt + 0.45, -5, axis="axis cs"),
			],
			xlabel = L"\lambda",
			ylabel = L"\mathcal{D}(\lambda)",
			xmin=0, xmax=5,
			ymin=-7, ymax=0,
			width = "9cm",
		)
	)
end

plot(p)

# +
using Optim
import LinearAlgebra: norm

p = let

	G = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left,
	                              xticklabels at=edge bottom, yticklabels at=edge left,
	                              horizontal sep=0.25cm, vertical sep=0.25cm",
	                  style="xlabel=\$x_1\$, ylabel=\$x_2\$")

	f = x -> norm(x) + sin(4atan(x[2], x[1]))
	g = x -> 2 - x[1]^2 - x[2]^2

	P = x -> max(g(x), 0)^2

	xdomain = (-3, 3)
	ydomain = (-3, 3)

	x_arr = Vector{Float64}[]
	x = [2.5, 2.5]
	push!(x_arr, x)

	for r in [0.0, 1.0, 2.0, 3.0]
	    obj = x -> f(x) + r*P(x)
	    x = optimize(obj, x+[0.25,0.25]).minimizer
	    push!(x_arr, x)

	    push!(G,
	        Axis([
	            Plots.Contour(obj, xdomain, ydomain, labels=false, xbins=151, ybins=151),
	            Plots.Contour(g, xdomain, ydomain, contour_style="draw color=black,labels=false", levels=[0.0]),
	            Plots.Scatter([x[1]], [x[2]], style="mark=*, mark size=1, mark options={draw=black, fill=black}"),
	            Plots.Node(L"x", x[1]+0.2, x[2]+0.2),
	            ], width="5.25cm", height="5.25cm", style="view={0}{90}", title="\$\\rho=$(round(Int, r))\$",
	        )
	    )
	end

	G
end

plot(p)

# +
using Optim

p = let

	f = x -> norm(x) + sin(4atan(x[2], x[1]))
	h = x -> 2 - x[1]^2 - x[2]^2

	P = x -> h(x)^2

	xs_P = Vector{Float64}[]
	x = x_start = [3, 2.5]
	c = c_start = 1.0
	γ = 1.5

	is = collect(1:10)

	for i in is
	    obj = x -> f(x) + c*P(x)
	    x = optimize(obj, x+randn(2)/5).minimizer
	    push!(xs_P, x)
	    c *= γ
	end

	hs = [h]
	c = c_start
	xs_ALM = Vector{Float64}[]
	λs = zeros(length(hs))
	x = x_start

	for i in is
	    P = x -> f(x) + c/2*sum(h(x)^2 for h in hs) -
	                        sum(λ*h(x) for (λ,h) in zip(λs, hs))
	    x = optimize(x -> f(x) + P(x), x+randn(2)/10).minimizer
	    λs -= c*[h(x) for h in hs]
	    c *= γ
	    push!(xs_ALM, x)
	end

	ax = Axis([
	        Plots.Linear(is, [h(x) for x in xs_P], style="black, solid, mark=none", legendentry="quadratic penalty method"),
	        Plots.Linear(is, [h(x) for x in xs_ALM], style="pastelBlue, solid, mark=none", legendentry="augmented Lagrangian"),
	      ], ymode="log", width="8cm", height="6cm", xlabel="number of iterations", ylabel=L"h(x)", style="legend cell align=left")
	ax.legendStyle = "{at={(1.05,1.0)},anchor=north west}"
	ax
end

plot(p)

# +
using Optim

p = let

	G = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left,
	                              xticklabels at=edge bottom, yticklabels at=edge left,
	                              horizontal sep=0.25cm, vertical sep=0.25cm",
	                  style="xlabel=\$x_1\$, ylabel=\$x_2\$")

	f = x -> norm(x) + sin(4atan(x[2], x[1]))
	g = x -> 2 - x[1]^2 - x[2]^2

	p = x -> g(x) <= 0 ? -1/g(x) : Inf

	xdomain = (-2.5, 2.5)
	ydomain = (-2.5, 2.5)

	x_arr = Vector{Float64}[]
	x = [3, 2.5]
	push!(x_arr, x)

	for ρ in [0.5, 1.0, 10.0, 100.0]
	    obj = x -> f(x) + p(x)/ρ
	    x = optimize(obj, x).minimizer
	    push!(x_arr, x)

	    push!(G,
	        Axis([
	            Plots.Contour(x->min(obj(x), 10.0), xdomain, ydomain, labels=false),
	            Plots.Contour(g, xdomain, ydomain, contour_style="draw color=black,labels=false", levels=[0.0]),
	            Plots.Scatter([x[1]], [x[2]], style="mark=*, mark size=1, mark options={draw=black, fill=black}"),
	            Plots.Node(L"x", x[1]+0.2, x[2]+0.2),
	            ], width="5.25cm", height="5.25cm", style="view={0}{90}", title="\$\\rho=$ρ\$",
	        )
	    )
	end

	G
end

plot(p)

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

# # Math Concepts
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +

p = let
	f = x -> cos(x) + x + 10x^1.5 + 3x^2
	f0 = x -> 110abs(cos(x))
	f1 = x -> 16x^1.5
	f2 = x -> 12x^2
	f3 = x -> 6x^3

	g = GroupPlot(2, 1, groupStyle="horizontal sep=1cm")
	a, b = 0, 2
	push!(g,
		Axis([
	    	Plots.Linear(f, (a, b), style="solid, thick, black, mark=none"),
	    	# Plots.Linear(f0, (a, b), style="solid, mark=none, pastelGreen"),
        	Plots.Linear(f1, (a, b), style="solid, mark=none, pastelBlue"),
        	Plots.Linear(f2, (a, b), style="solid, mark=none, pastelPurple"),
        	Plots.Linear(f3, (a, b), style="solid, mark=none, pastelSeaGreen"),
        ], xmin=a, xmax=b, xlabel=L"x", width="6cm")
	)

	a,b = 1,1000
	push!(g,
		Axis([
	    	Plots.Linear(f, (a, b),  style="solid, thick, black, mark=none", legendentry=L"f(x)"),
        	# Plots.Linear(f0, (a, b), style="solid, mark=none, pastelGreen", legendentry=L"110|\cos(x)|"),
        	Plots.Linear(f1, (a, b), style="solid, mark=none, pastelBlue", legendentry=L"16x^{3/2}"),
        	Plots.Linear(f2, (a, b), style="solid, mark=none, pastelPurple", legendentry=L"12x^2"),
        	Plots.Linear(f3, (a, b), style="solid, mark=none, pastelSeaGreen", legendentry=L"6x^3"),
        ], xmin=a, xmax=b, xlabel=L"x", width="6cm", style="xmode=log, ymode=log, legend style={at={(-0.11,-0.27)},anchor=north,legend columns=5}")
	)
	g
end

plot(p)

# +

p = let
	f = x -> cos(x)
	f′ = x -> -sin(x)
	f′′ = x ->-cos(x)
	f′′′ = x ->sin(x)
	f′′′′ = x ->cos(x)
	f′′′′′ = x ->-sin(x)

	p = 1.0
	a = p-5.0
	b = p+5.0

	taylor0 = x -> f(p)
	taylor1 = x -> f(p) + f′(p)*(x-p)
	taylor2 = x -> f(p) + f′(p)*(x-p) + f′′(p)*(x-p)^2/2
	taylor3 = x -> f(p) + f′(p)*(x-p) + f′′(p)*(x-p)^2/2 + f′′′(p)*(x-p)^3/6
	taylor4 = x -> f(p) + f′(p)*(x-p) + f′′(p)*(x-p)^2/2 + f′′′(p)*(x-p)^3/6 + f′′′′(p)*(x-p)^4/24
	taylor5 = x -> f(p) + f′(p)*(x-p) + f′′(p)*(x-p)^2/2 + f′′′(p)*(x-p)^3/6 + f′′′′(p)*(x-p)^4/24 + f′′′′′(p)*(x-p)^5/120

	x_arr = collect(range(a, stop=b, length=151))

	Axis([
	    Plots.Linear(x_arr, f.(x_arr), style="solid, thick, black, mark=none", legendentry=L"\cos(x)"),
        Plots.Linear(x_arr, taylor0.(x_arr), style="solid, mark=none", legendentry="0th degree"),
        Plots.Linear(x_arr, taylor1.(x_arr), style="solid, mark=none", legendentry="1st degree"),
        Plots.Linear(x_arr, taylor2.(x_arr), style="solid, mark=none", legendentry="2nd degree"),
        Plots.Linear(x_arr, taylor3.(x_arr), style="solid, mark=none", legendentry="3rd degree"),
        Plots.Linear(x_arr, taylor4.(x_arr), style="solid, mark=none", legendentry="4th degree"),
        Plots.Linear(x_arr, taylor5.(x_arr), style="solid, mark=none", legendentry="5th degree"),
        Plots.Scatter([p], [f(p)], style="solid, mark=*, mark size=1, mark options={draw=none, fill=black}"),
        ], xmin=a, xmax=b, ymin=-3, ymax=3, xlabel=L"x", width="10cm", height="8cm", style="legend columns=1, legend pos=outer north east, cycle list name = pastelcolors")
end

plot(p)

# +
using Distributions

p = let

	N = Normal(0.0,1.0)

	Axis(Plots.Linear(x->pdf(N, x), (-3,3), style="black,solid,mark=none"),
		 xmin=-3, xmax=3, ymin=0, ymax=0.6, width="9cm",
		 xlabel=L"x", ylabel=L"p(x\mid \mu, \sigma)",
		 style="xtick={-1,0,1}, ytick=\\empty, xticklabels={\$\\mu-\\sigma\$,\$\\mu\$,\$\\mu+\\sigma\$}, xticklabel style={text height=2ex}",
		 )
end

plot(p)

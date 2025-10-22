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

# # Test Functions
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
p = let
	function ackley(x, a=20, b=0.2, c=2π)
		d = length(x)
		return -a*exp(-b*sqrt(sum(x.^2)/d)) -
		          exp(sum(cos.(c*xi) for xi in x)/d) + a + ℯ
	end

	xdomain = (-30, 30)
	ydomain = (-30, 30)

	plots = Plots.Plot[]
	push!(plots, Plots.Image((x,y)->ackley([x,y]), xdomain, ydomain, xbins=600, ybins=600, colormap = pasteljet, colorbar = false))
	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="axis equal, view={0}{90}")
end

plot(p)

# +
p = let
	booth(x) = (x[1]+2x[2]-7)^2 + (2x[1]+x[2]-5)^2

	xdomain = (-10, 10)
	ydomain = (-10, 10)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(booth, xdomain, ydomain, style="width=\\columnwidth", xbins=100, ybins=100, levels=[0.5,2,5,10,20,50,100,500,1000]))
	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +

p = let
	branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = a*(x[2] - b*x[1]^2 + c*x[1] - r)^2 + s*(1-t)*cos(x[1]) + s

	xdomain = (-2π, 6π)
	ydomain = ( -π, 7π)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(branin, xdomain, ydomain, levels=[1,2,3,5,10,20,50,100], style="width=\\columnwidth", xbins=100, ybins=100))
	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +
using LinearAlgebra
p = let
	kochenderfer(x; a=1, b=1, c=4) = a*norm(x) + b*sin(c*atan(x[2], x[1]))

	xdomain = ( -3, 3)
	ydomain = ( -3, 3)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(kochenderfer, xdomain, ydomain, levels=collect(-4:4), style="width=\\columnwidth", xbins=101, ybins=101))
	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +
p = let
michalewicz(x; m=10) = -sum(sin(v)*sin(i*v^2/π)^(2m) for (i,v) in enumerate(x))

xdomain = (0, 4)
ydomain = (0, 4)

plots = Plots.Plot[]
# push!(plots, Plots.Contour(michalewicz, xdomain, ydomain, levels=[-1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2], style="width=\\columnwidth", xbins=101, ybins=101))
# Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
Axis(Plots.Image((x,y)->michalewicz([x,y]), xdomain, ydomain, xbins=600, ybins=600, colormap = pasteljet, colorbar = false),
xmin=xdomain[1], xmax=xdomain[2], ymin=ydomain[1], ymax=ydomain[2], width="8cm", height="8cm", style="view={0}{90}", xlabel=L"x_1", ylabel=L"x_2")
end

plot(p)

# +
p = let
	rosenbrock(x; a=1, b=5) = (a-x[1])^2 + b*(x[2] - x[1]^2)^2

	xdomain = (-2, 2)
	ydomain = (-2, 2)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(rosenbrock, xdomain, ydomain, levels=[1,2,3,5,10,20,50,100], style="width=\\columnwidth", xbins=100, ybins=100))
	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +
p = let
	wheeler(x, a=1.5) = -exp(-(x[1]*x[2] - a)^2 -(x[2]-a)^2)

	xdom = (-8,25)
	ydom = (-3,6)

	Axis(Plots.Image((x,y)->wheeler([x,y]), xdom, ydom, xbins=1200, ybins=480, colormap = pasteljet, colorbar = false),
	     xmin=xdom[1], xmax=xdom[2], ymin=ydom[1], ymax=ydom[2], width="12cm", height="4.8cm", style="view={0}{90}", xlabel=L"x_1", ylabel=L"x_2")
end

plot(p)

# +
p = let

	f = x -> -exp(-(x[1]*x[2] - 1.5)^2 -(x[2]-1.5)^2)

	xdomain = (0, 3)
	ydomain = (0, 3)

	p_contour = Plots.Contour(f, xdomain, ydomain, style="width=\\columnwidth", xbins=151, ybins=151) # , levels=[-15,-10,-5,-1,0,0.5,1,1.5,2,2.5,3]

	plots = Plots.Plot[]
	push!(plots, p_contour)
	ax = Axis(plots, ymin=0, xmin=0, ymax=3, xmax=3, width="9cm", height="9cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

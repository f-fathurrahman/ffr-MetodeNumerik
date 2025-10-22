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

# # First Order
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Vec

p = let
	function secant_method(df, x1, x2, ϵ)
	    df1 = df(x1)

	    delta = Inf
	    while abs(delta) > ϵ
	    	df2 = df(x2)
	        delta = (x2 - x1)/(df2 - df1)*df2
	        x1, x2, df1 = x2, x2 - delta, df2
	    end
	    x2
	end

	f = (x,y) -> x^2 + y^2 + y*sin(x)
	df = (x,y) -> [2x + y*cos(x), 2y + sin(x)]

	xdomain = (-2, 4)
	ydomain = (-2, 4)

	p0 = (-0.5, 3.5)
	pts = Tuple{Float64,Float64}[p0]
	for i in 1 : 3
		x,y = pts[end]
		dp = df(x, y)
		vdp = VecE2{Float64}(-dp[1], -dp[2])

		f1d = a -> begin
			x2 = x + a*dp[1]
			y2 = y + a*dp[2]

			da = df(x2, y2)
			pa = VecE2{Float64}(da[1], da[2])
			proj(pa, vdp, Float64)
		end
		alpha = secant_method(f1d, 0.0, 1.0, 0.001)

		push!(pts, (x + alpha*dp[1], y + alpha*dp[2]))
	end

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=0:2:10, style="width=\\columnwidth"))
	push!(plots, Plots.Linear3([p[1] for p in pts], [p[2] for p in pts], [f(p[1], p[2]) for p in pts], style="black, solid, mark=*, mark size = 1, mark options={draw=black}"))

	g = GroupPlot(2, 1)
	push!(g, Axis(plots, height="5cm", xlabel=L"x_1", ylabel=L"x_2", zlabel=L"y", style="xtick=\\empty, ytick=\\empty, contour/labels=false"))
	push!(g, Axis(plots, width="5cm", height="5cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}"))
	g
end

plot(p)

# +
using Vec
using LinearAlgebra

pts_gradient_banana = Tuple{Float64,Float64}[]

p = let
	function secant_method(df, x1, x2, ϵ)
	    df1 = df(x1)

	    delta = Inf
	    while abs(delta) > ϵ
	    	df2 = df(x2)
	        delta = (x2 - x1)/(df2 - df1)*df2
	        x1, x2, df1 = x2, x2 - delta, df2
	    end
	    x2
	end


	f = (x,y) -> (1-x)^2 + 5*(y - x^2)^2
	df = (x,y) -> [2*(10*x^3-10*x*y+x-1), 10*(y-x^2)]

	xdomain = (-2, 2)
	ydomain = (-2, 2)

	p0 = (-1, -1)
	pts = Tuple{Float64,Float64}[p0]
	for i in 1 : 10
		x,y = pts[end]
		dp = normalize(-VecE2{Float64}(df(x, y)...))

		f1d = a -> begin
			x2 = x + a*dp.x
			y2 = y + a*dp.y

			da = df(x2, y2)
			pa = VecE2{Float64}(da[1], da[2])
			proj(pa, dp, Float64)
		end
		alpha = secant_method(f1d, 0.0, 1.0, 0.0001)
		push!(pts, (x + alpha*dp.x, y + alpha*dp.y))
	end

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[1,2,3,5,10,20,50,100], style="width=\\columnwidth"))
	push!(plots, Plots.Linear3([p[1] for p in pts], [p[2] for p in pts], [f(p[1], p[2]) for p in pts], style="black, solid, mark=none"))

	append!(pts_gradient_banana, pts)

	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)
# -

abstract type DescentMethod end

# +
p = let

	A = Float64[1 -0.9; -0.9 1]
	b = Float64[0, 0]
	f = x -> (x'*A*x)[1]
	∇f = x -> 2*A'*x

	x = Float64[-1.1, -1.8]
	pts = Vector{Float64}[x]

	function get_optimal_step_size(d, A, b)
		return -(d'*(A*x+b))[1]/(d'*A*d)[1]
	end

	# FIRST ITERATION
	d = -∇f(x)
	α = get_optimal_step_size(d, A, b)
	x += α*d
	push!(pts, x)

	# SECOND ITERATION
	g = ∇f(x)
	β = (g'*A*d)[1] / (d'*A*d)[1]
	d = -g + β*d
	α = get_optimal_step_size(d, A, b)
	x += α*d
	push!(pts, x)

	xdomain = (-2, 2)
	ydomain = (-2, 2)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain))
	push!(plots, Plots.Linear3([x[1] for x in pts],
	                           [x[2] for x in pts],
	                           [f(x) for x in pts], style="black, solid, mark=*, mark size = 1, mark options={draw=black}"))
	push!(plots, Plots.Node(L"1", pts[1][1], pts[1][2], style="right"))
	push!(plots, Plots.Node(L"2", pts[2][1], pts[2][2], style="left"))
	push!(plots, Plots.Node(L"3", pts[3][1], pts[3][2], style="right"))


	Axis(plots, width="\\columnwidth", height="\\columnwidth", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)

# +
using Vec

p = let
	function secant_method(df, x1, x2, ϵ)
	    df1 = df(x1)

	    delta = Inf
	    while abs(delta) > ϵ
	    	df2 = df(x2)
	        delta = (x2 - x1)/(df2 - df1)*df2
	        x1, x2, df1 = x2, x2 - delta, df2
	    end
	    x2
	end

	f = x -> (1-x[1])^2 + 5*(x[2] - x[1]^2)^2
	df = x -> [2*(10*x[1]^3-10*x[1]*x[2]+x[1]-1), 10*(x[2]-x[1]^2)]

	x0 = Float64[-1, -1]
	pts = Vector{Float64}[x0]
	g0 = df(x0)
	d0 = -g0
	f1d = a -> begin
		da = df(x0 + a*d0)
		proj(VecE2{Float64}(da[1], da[2]), VecE2{Float64}(d0[1], d0[2]), Float64)
	end
	alpha = secant_method(f1d, 0.0, 0.1, 0.0001) # NOTE: is sensitive to the value of 0.1, does not work with 1.0!
	push!(pts, x0 + alpha*d0)

	for i in 1 : 5
		x = pts[end]
		g1 = df(x)
		β = max(0.0, dot(g1, g1-g0)/dot(g0,g0))
		d1 = -g1 + β*d0

		f1d = a -> begin
			da = df(x + a*d1)
			proj(VecE2{Float64}(da[1], da[2]), VecE2{Float64}(d1[1], d1[2]), Float64)
		end
		alpha = secant_method(f1d, 0.0, 0.1, 0.0001)
		push!(pts, x + alpha*d1)

		d0 = d1
		g0 = g1
	end

	xdomain = (-2, 2)
	ydomain = (-2, 2)

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[1,2,3,5,10,20,50,100], style="width=\\columnwidth"))
	push!(plots, Plots.Linear3([p[1] for p in pts_gradient_banana],
	                           [p[2] for p in pts_gradient_banana],
	                           [f([p[1], p[2]]) for p in pts_gradient_banana], style="gray, solid, mark=none"))
	push!(plots, Plots.Linear3([p[1] for p in pts], [p[2] for p in pts], [f(p) for p in pts], style="black, solid, mark=none"))


	Axis(plots, width="8cm", height="8cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)
# -

p = let
	Axis([
		Plots.Linear(x->-exp(-x^2), (-4,4), style="solid, black, mark=none"),
		Plots.Linear(x->-exp(-x^2), (-4,-2.5), style="solid, thick, pastelBlue, mark=none"),
		Plots.Linear(x->-exp(-x^2), (2.5,4), style="solid, thick, pastelBlue, mark=none"),
		], xlabel=L"x", ylabel=L"-\exp(-x^2)", width="9cm", style="enlarge x limits=0"
	)
end
plot(p)

# +
using Vec

struct GradientDescent <: DescentMethod
	α
end
init!(M::GradientDescent, f, ∇f, x) = M
function step!(M::GradientDescent, f, ∇f, x)
	α, g = M.α, ∇f(x)
	return x - α*g
end
mutable struct Momentum <: DescentMethod
	α # learning rate
	β # momentum decay
	v # momentum
end
function init!(M::Momentum, f, ∇f, x)
	M.v = zeros(length(x))
	return M
end
function step!(M::Momentum, f, ∇f, x)
	α, β, v, g = M.α, M.β, M.v, ∇f(x)
	v[:] = β*v - α*g
	return x + v
end

p = let
	f = x -> (1-x[1])^2 + 100*(4x[2] - x[1]^2)^2
	∇f = x -> [2*(200x[1]^3 - 800x[1]*x[2] + x[1] - 1), -800*(x[1]^2 - 4x[2])]

	xdomain = (-3, 2)
	ydomain = (-0.5, 2)

	function this_step!(M::DescentMethod, v::VecE2{Float64})
	    x = Float64[v.x, v.y]
	    return VecE2{Float64}(step!(M, f, ∇f, x)...)
	end
	function run_descent_method(M::DescentMethod, x₀::VecE2{Float64}, N::Int)
	    pts = [x₀]
	    init!(M, f, ∇f, Float64[x₀.x, x₀.y])
	    for i in 1 : N
	        push!(pts, this_step!(M, pts[end]))
	    end
	    return pts
	end
	function get_descent_plot(pts::Vector{VecE2{Float64}}, name, color::String="black")
	    Plots.Linear([P.x for P in pts], [P.y for P in pts], style="thick, $color, solid, mark=none, line join=round", legendentry=name) #$
	end

	x₀ = VecE2{Float64}(-2,1.5)
	N = 40

	stuff = Tuple{DescentMethod, String, String}[]
	push!(stuff, (GradientDescent(0.0003), "gradient descent", "pastelRed"))
	push!(stuff, (Momentum(0.0003, 0.9, zeros(2)), "momentum", "black"))

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[2,10,50,200,500], style="forget plot, width=\\columnwidth", xbins=150, ybins=150))
	for (M, name, color) in stuff
	    pts = run_descent_method(M, x₀, N)
	    push!(plots, get_descent_plot(pts, name, color))
	end
	Axis(plots, width="12cm", height="6cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, view={0}{90}, legend pos=outer north east")
end
plot(p)

# +
using Vec

mutable struct NesterovMomentum <: DescentMethod
	α # learning rate
	β # momentum decay
	v # momentum
end
function init!(M::NesterovMomentum, f, ∇f, x)
	M.v = zeros(length(x))
	return M
end
function step!(M::NesterovMomentum, f, ∇f, x)
	α, β, v = M.α, M.β, M.v
	v[:] = β*v - α*∇f(x + β*v)
	return x + v
end

p = let
	f = x -> (1-x[1])^2 + 100*(4x[2] - x[1]^2)^2
	∇f = x -> [2*(200x[1]^3 - 800x[1]*x[2] + x[1] - 1), -800*(x[1]^2 - 4x[2])]

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

	x₀ = VecE2{Float64}(-2,1.5)
	dm = GradientDescent(0.0003)
	N = 40

	stuff = Tuple{DescentMethod, String, String}[]
	# push!(stuff, (GradientDescent(0.0003), "gradient descent", "black"))
	push!(stuff, (Momentum(0.0003, 0.9, zeros(2)), "momentum", "black"))
	push!(stuff, (NesterovMomentum(0.0002, 0.92, zeros(2)), "Nesterov momentum", "pastelRed"))
	# push!(stuff, (Adagrad(0.1, 1e-8, zeros(2)), "Adagrad", "green"))
	# push!(stuff, (RMSProp(0.65, 0.45, 1e-8, zeros(2)), "RMSProp", "green"))
	# push!(stuff, (Adadelta(0.8, 0.8, 1e-3, zeros(2), zeros(2)), "Adadelta", "orange"))
	# push!(stuff, (Adam(0.8, 0.9, 0.9, 1e-8, 0, zeros(2), zeros(2)), "Adam", "purple"))

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[2,10,50,200,500], style="forget plot, width=\\columnwidth", xbins=150, ybins=150))
	for (M, name, color) in stuff
	    pts = run_descent_method(M, x₀, N)
	    push!(plots, get_descent_plot(pts, name, color))
	end
	Axis(plots, width="12cm", height="6cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, view={0}{90}, legend pos=outer north east")
end
plot(p)

# +
using Vec

mutable struct HyperGradientDescent <: DescentMethod
	α0 # initial learning rate
	μ # learning rate of the learning rate
	α # current learning rate
	g_prev # previous gradient
end
function init!(M::HyperGradientDescent, f, ∇f, x)
	M.α = M.α0
	M.g_prev = zeros(length(x))
	return M
end
function step!(M::HyperGradientDescent, f, ∇f, x)
	α, μ, g, g_prev = M.α, M.μ, ∇f(x), M.g_prev
	α = α + μ*(g⋅g_prev)
	M.g_prev, M.α = g, α
	return x - α*g
end

mutable struct HyperNesterovMomentum <: DescentMethod
	α0 # initial learning rate
	μ # learning rate of the learning rate
	β # momentum decay
	v # momentum
	α # current learning rate
	g_prev # previous gradient
end
function init!(M::HyperNesterovMomentum, f, ∇f, x)
	M.α = M.α0
	M.v = zeros(length(x))
	M.g_prev = zeros(length(x))
	return M
end
function step!(M::HyperNesterovMomentum, f, ∇f, x)
	α, β, μ = M.α, M.β, M.μ
	v, g, g_prev = M.v, ∇f(x), M.g_prev
	α = α - μ*(g⋅(-g_prev - β*v))
	v[:] = β*v + g
	M.g_prev, M.α = g, α
	return x - α*(g + β*v)
end

p = let
	f = x -> (1-x[1])^2 + 100*(4x[2] - x[1]^2)^2
	∇f = x -> [2*(200x[1]^3 - 800x[1]*x[2] + x[1] - 1), -800*(x[1]^2 - 4x[2])]

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

	x₀ = VecE2{Float64}(-2,1.5)
	dm = GradientDescent(0.0003)
	N = 40

	stuff = Tuple{DescentMethod, String, String}[]
	push!(stuff, (HyperGradientDescent(0.0004, 8e-13, NaN, zeros(2)), "hypergradient", "black"))
	push!(stuff, (HyperNesterovMomentum(0.00023, 1e-12, 0.93, zeros(2), NaN, zeros(2)), "hyper-Nesterov", "pastelRed"))

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[2,10,50,200,500], style="forget plot, width=\\columnwidth", xbins=150, ybins=150))
	for (M, name, color) in stuff
	    pts = run_descent_method(M, x₀, N)
	    push!(plots, get_descent_plot(pts, name, color))
	end
	Axis(plots, width="12cm", height="6cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, view={0}{90}, legend pos=outer north east")
end
plot(p)

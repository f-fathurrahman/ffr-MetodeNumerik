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

# # Surrogate Models
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Random
using LinearAlgebra

function design_matrix(X)
	n, m = length(X[1]), length(X)
	return [j==0 ? 1. : X[i][j] for i in 1:m, j in 0:n]
end
function linear_regression(X, y)
	θ = pinv(design_matrix(X))*y
	return x -> θ⋅[1; x]
end

p = let

	G = GroupPlot(2,2,groupStyle="horizontal sep=01cm, vertical sep=1.5cm", style="xlabel=\$x\$, ylabel=\$y\$, width=6cm, ytick=\\empty, xtick=\\empty, enlargelimits=0, ymin=-0.2, ymax=0.4, xmin=0, xmax=1")

	Random.seed!(0)
	f = x -> 0.2x + 0.2randn()

	function get_plot(X, ys, title)
		fhat = linear_regression(X, ys)
		plots = Plots.Plot[]
		push!(plots, Plots.Scatter(X, ys, style="mark=*, mark size=1, mark options={draw=black, fill=black}"))
		for (x, y)  in zip(X, ys)
			dist = abs(y - fhat(x))
			if dist > 0.01
				push!(plots, Plots.Linear([x,x], [y, fhat(x)], style="solid, gray, mark=none"))
			end
		end
		push!(plots, Plots.Linear([0,1], [fhat(0), fhat(1)], style="solid, pastelBlue, mark=none"))
		Axis(plots, title=title)
	end

	X = [0.2]
	ys = f.(X)
	push!(G, get_plot(X, ys, L"m < n+1"))

	push!(X, 0.6)
	push!(ys, f(X[end]))
	push!(G, get_plot(X, ys, L"m = n+1"))

	push!(G, get_plot([X[1], X[1]], [ys[1], f(X[1])], "nonindependent points"))

	append!(X, [0.1, 0.35, 0.4, 0.4, 0.8])
	append!(ys, f.(X[end-4:end]))
	push!(G, get_plot(X, ys, L"m > n+1"))

	G
end

plot(p)

# +
using Random
using LinearAlgebra

polynomial_bases_1d(i, k) = [x->x[i]^p for p in 0:k]
function polynomial_bases(n, k)
	bases = [polynomial_bases_1d(i, k) for i in 1 : n]
	terms = Function[]
	for ks in Iterators.product([0:k for i in 1:n]...)
		if sum(ks) ≤ k
			push!(terms,
				x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
		end
	end
	return terms
end
function regression(X, y, bases)
    B = [b(x) for x in X, b in bases]
    θ = pinv(B)*y
    return x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
end

p = let

	a = 45
	b = 180+123
	G = GroupPlot(2,1, groupStyle="horizontal sep = 1cm", style="width=6cm")

	f = x -> x^2 + 0.1
	Random.seed!(0)
	X = [0.1, 0.2, 0.35, 0.4, 0.6, 0.8]
	y = f.(X)

	fhat = regression(X, y, polynomial_bases(1,2))

	push!(G, Axis([
	        Plots.Scatter(X, y, style="mark size=1, mark options={draw=black, fill=black}")
	        Plots.Linear(x->fhat([x]), (0,1), style="solid, pastelBlue, mark=none")
	        ], ylabel=L"y", xlabel=L"x", style="xticklabels=\\empty, yticklabels=\\empty"))

	x_arr = collect(range(0, stop=1, length=101))
	x2_arr = x_arr.^2
	y_arr = [fhat(x) for x in x_arr]

	push!(G, Axis([
	        Plots.Command("\\addplot3[patch,patch type=rectangle, faceted color=none, color=pastelBlue!40] coordinates {(0,0,0.1) (1,0,0.1) (1,1,1.1) (0,1,1.1)}"),
	        Plots.Linear3(X, X.^2, y, mark="*", style="mark size=1.1, only marks, mark options={draw=black, fill=black}", legendentry=L"f(x)"),
	        Plots.Linear3(x_arr, x2_arr, y_arr, style="pastelBlue, solid, mark=none", legendentry=L"\hat{f}(x)"),
	    ], xlabel=L"x", ylabel=L"x^2", zlabel=L"y", width="6cm", style="view={45}{20}, axis equal, xticklabels=\\empty, yticklabels=\\empty, zticklabels=\\empty", legendPos="outer north east"))

	G
end

plot(p)

# +
using Random

function sinusoidal_bases_1d(j, k, a, b)
	T = b[j] - a[j]
	bases = Function[x->1.0]
	for i in 1 : k
		push!(bases, x->sin(2π*i*x[j]/T))
		push!(bases, x->cos(2π*i*x[j]/T))
	end
	return bases
end
function sinusoidal_bases(k, a, b)
	n = length(a)
	bases = [sinusoidal_bases_1d(i, k, a, b) for i in 1 : n]
	terms = Function[]
	for ks in Iterators.product([0:2k for i in 1:n]...)
		powers = [div(k+1,2) for k in ks]
		if sum(powers) ≤ k
			push!(terms,
				x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
		end
	end
	return terms
end

p = let

	Random.seed!(0)
	a = [0.0]
	b = [10.0]
	f = x -> sin(x) + 0.2randn()

	G = GroupPlot(2,2,groupStyle="horizontal sep=01cm, vertical sep=1.5cm", style="xlabel=\$x\$, ylabel=\$y\$, width=6cm, ytick=\\empty, xtick=\\empty, enlargelimits=0, ymin=-2, ymax=2, xmin=$(a[1]), xmax=$(b[1])")

	function get_plot(X, ys, title)
		fhat = regression(X, ys, sinusoidal_bases(1, a, b))
		plots = Plots.Plot[]
		push!(plots, Plots.Scatter(X, ys, style="mark=*, mark size=1, mark options={draw=black, fill=black}"))
		for (x, y)  in zip(X, ys)
			dist = abs(y - fhat(x))
			if dist > 0.01
				push!(plots, Plots.Linear([x,x], [y, fhat([x])], style="solid, gray, mark=none"))
			end
		end
		push!(plots, Plots.Linear(x->fhat(x), (a[1],b[1]), style="solid, pastelBlue, mark=none"))
		Axis(plots, title=title)
	end

	X = [2.0]
	ys = f.(X)
	push!(G, get_plot(X, ys, L"m < n+1"))

	push!(X, 5.6)
	push!(ys, f(X[end]))
	push!(G, get_plot(X, ys, L"m = n+1"))

	push!(G, get_plot([X[1], X[1]], [ys[1], f(X[1])], "nonindependent points"))

	append!(X, [1.1, 9.35, 4, 4.4, 8.8])
	append!(ys, f.(X[end-4:end]))
	push!(G, get_plot(X, ys, L"m > n+1"))

	G
end

plot(p)

# +
p = let
	xdom = (0,2)
	ydom = (-2,2)

	g = GroupPlot(3,2, groupStyle="horizontal sep=1.75cm, vertical sep=1cm, xlabels at=edge bottom, ylabels at=edge left, xticklabels at=edge bottom, yticklabels at=edge bottom", style="xlabel=\$r\$, ylabel=\$\\psi\$, title style={text height=2ex}")

	function add_plot!(f, title)
	    push!(g,
	        Axis(
	            Plots.Linear(r -> f(abs(r)), xdom, style="thick, solid, black, mark=none"),
	            width="5cm", height="5cm", title=title,
	        ),
	    )
	end

	add_plot!(r -> r, "linear: \$r\$")
	add_plot!(r -> r^3, "cubic: \$r^3\$")
	add_plot!(r -> r^2 * log(r), "thin plate spline: \$r^2 \\log r\$")
	add_plot!(r -> exp(-r^2), "Gaussian: \$e^{-r^2/2\\sigma^2} \$")
	add_plot!(r -> (r^2 + 1)^0.5, "multiquadratic: \$(r^2 + \\sigma^2)^\\frac{1}{2}\$")
	add_plot!(r -> (r^2 + 1)^-0.5, "inverse multiquadratic: \$(r^2 + \\sigma^2)^{-\\frac{1}{2}}\$")

	g
end

plot(p)

# +
using Random

radial_bases(ψ, C, p=2) = [x->ψ(norm(x - c, p)) for c in C]

p = let

	f = x -> x * sin(5x)

	Random.seed!(6)
	X = [rand() for i in 1 : 4]
	y = (x->f(x[1])).(X)

	Axis([
	        Plots.Linear(f, (0,1), style="thick, solid, black, mark=none", legendentry=L"x \sin(5x)"),
	        Plots.Linear(regression(X, y, radial_bases(r -> exp(-2*r^2), X)), (0,1), style="thick, pastelBlue!40", legendentry=L"\psi = \exp(-2r^2)"),
	        Plots.Linear(regression(X, y, radial_bases(r -> exp(-5*r^2), X)), (0,1), style="thick, pastelBlue!70", legendentry=L"\psi = \exp(-5r^2)"),
	        Plots.Linear(regression(X, y, radial_bases(r -> exp(-10*r^2), X)), (0,1), style="thick, pastelBlue", legendentry=L"\psi = \exp(-10r^2)"),
	        Plots.Scatter(X, y, style="black, mark=*, mark size=1, mark options={draw=none, fill=black}"),
	        ], xlabel=L"x", ylabel=L"y", style="enlarge x limits=0", width="8cm", legendPos="outer north east"
	    )
end

plot(p)

# +
using Random

function regression(X, y, bases, λ)
    B = [b(x) for x in X, b in bases]
    θ = (B'B + λ^2*I)\B'y
    return x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
end

p = let

	f = x -> x * sin(5x)

	Random.seed!(6)
	X = rand(10)
	y = f.(X) + randn(length(X))/10

	Axis([
	        Plots.Linear(f, (0,1), style="thick, solid, black, mark=none", legendentry=L"x \sin (5x)"),
	        Plots.Linear(regression(X, y, radial_bases(r -> exp(-5*r^2), X)), (0,1), style="thick, pastelBlue", legendentry=L"\lambda = 0"),
	        Plots.Linear(regression(X, y, radial_bases(r -> exp(-5*r^2), X), 0.1), (0,1), style="thick, pastelPurple", legendentry=L"\lambda = 0.1"),
	        Plots.Linear(regression(X, y, radial_bases(r -> exp(-5*r^2), X), 0.5), (0,1), style="thick, pastelRed", legendentry=L"\lambda = 0.5"),
	        Plots.Scatter(X, y, style="black, mark=*, mark size=1, mark options={draw=none, fill=black}"),
	        ], xlabel=L"x", ylabel=L"y", style="enlarge x limits=0", width="8cm", legendPos="outer north east", ymin=-2, ymax=2,
	    )
end

plot(p)

# +
p = let
	variance = x->log(1+exp(5(x-1))) + 0.1
	bias = x->variance(1-x)
	toterr = x->bias(x) + variance(x)

	Axis([
	        Plots.Linear(variance, (0,1), style="solid, pastelBlue, mark=none", legendentry="variance"),
	        Plots.Linear(bias, (0,1), style="solid, pastelRed, mark=none", legendentry=L"\text{bias}^2"),
	        Plots.Linear(toterr, (0,1), style="solid, black, mark=none", legendentry="total error"),
	        Plots.Linear([0.5,0.5],[0.1,0.9], style="solid, gray, mark=none", legendentry="optimal tradeoff"),
	    ],
	    width="9cm", xlabel="complexity", ylabel="error",
	    ymin=0, ymax=1, xmin=0, xmax=1,
	    style="axis lines=middle, xtick=\\empty, ytick=\\empty, legend cell align=left, legend style={draw=none, fill=none, at={(0.5,-0.25)}, anchor=north, legend columns=1}, x label style={at={(axis description cs:0.75,-0.05)},anchor=north},",
	)
end

plot(p)

# +
using Random
using Distributions
using LinearAlgebra

function regression(X, y, bases)
    B = [b(x) for x in X, b in bases]
    θ = pinv(B)*y
    return x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
end

f = x->2/(1+exp(-x))
r = 4

str = """
\\begin{tikzpicture}[x=0.5cm, y=0.5cm]
    \\coordinate (A) at (2,6);")
    \\coordinate (B) at (\$(A) + (4.5,0)\$);
    \\coordinate (C) at (\$(A) + (0,-4.5)\$);
    \\coordinate (D) at (C -| B);

    \\draw (A) circle (1.75);
    \\draw[gray] (\$(A) + (0:2)\$) -- ++(180:4);
    \\draw[gray] (\$(A) + (90:2)\$) -- ++(270:4);
    \\draw (B) circle (1.75);
    \\draw[gray] (\$(B) + (0:2)\$) -- ++(180:4);
    \\draw[gray] (\$(B) + (90:2)\$) -- ++(270:4);
    \\draw (C) circle (1.75);
    \\draw[gray] (\$(C) + (0:2)\$) -- ++(180:4);
    \\draw[gray] (\$(C) + (90:2)\$) -- ++(270:4);
    \\draw (D) circle (1.75);
    \\draw[gray] (\$(D) + (0:2)\$) -- ++(180:4);
    \\draw[gray] (\$(D) + (90:2)\$) -- ++(270:4);

    \\node[rotate=90, anchor=south] at (\$(A) + (-2.25,0)\$) {low bias\\vphantom{g}};
    \\node[rotate=90, anchor=south] at (\$(C) + (-2.25,0)\$) {high bias};
    \\node[anchor=south] at (\$(A) + (0,2.25)\$) {low variance\\vphantom{g}};
    \\node[anchor=south] at (\$(B) + (0,2.25)\$) {high variance};
    """

function print_coordinate(str, p, f1, f2, x1, x2)
    x = f1(x1) - f(x1)
    y = f2(x2) - f(x2)
    if hypot(x,y) < 3
        str *= "\n\\fill[red, opacity=0.5] (\$($p) + ($x, $y)\$) circle (0.05);"
    end
    return str
end

Random.seed!(0)
for i in 1 : 100

    X1 = collect(range(-r, stop=r, length=10))
    X2 = collect(range(-r, stop=r, length=10))
    y1 = (x->f(x) + 0.25randn()).(X1)
    y2 = (x->f(x) + 0.25randn()).(X2)

    f11 = regression(X1, y1, [x->1.0, x->x, x->x^3])
    f12 = regression(X2, y2, [x->1.0, x->x, x->x^3])
    f21 = regression(X1, y1, [x->1.0, x->x, x->x^2, x->x^3, x->x^4, x->x^5, x->x^6, x->x^7, x->x^8, x->x^9])
    f22 = regression(X2, y2, [x->1.0, x->x, x->x^2, x->x^3, x->x^4, x->x^5, x->x^6, x->x^7, x->x^8, x->x^9])
    f31 = regression(X1, y1, [x->x])
    f32 = regression(X2, y2, [x->x])
    f41 = regression(X1, y1, [x->sin(x/4), x->sin(x/8), x->sin(x/16), x->sin(x/32), x->sin(x/64)])
    f42 = regression(X2, y2, [x->sin(x/4), x->sin(x/8), x->sin(x/16), x->sin(x/32), x->sin(x/64)])

    x1 = rand(Uniform(-4,4))
    x2 = rand(Uniform(-4,4))

    str = print_coordinate(str, "A", f11, f12, x1, x2)
    str = print_coordinate(str, "B", f21, f22, x1, x2)
    str = print_coordinate(str, "C", f31, f32, x1, x2)
    str = print_coordinate(str, "D", f41, f42, x1, x2)
end

str *= "\\end{tikzpicture}"

using TikzPictures
pic = TikzPicture(str, preamble="\\usetikzlibrary{calc}")

# +
using Random

p = let
	g = GroupPlot(2,1)

	f = x -> x + x^2 + randn()/2
	Random.seed!(0)
	X = [0.1, 0.2, 0.35, 0.4, 0.6, 0.8]
	y = f.(X)

	g = GroupPlot(3,2)
	for k in 0 : 5
	    fhat = regression(X, y, polynomial_bases(1, k))
	    err = @sprintf("%.4f", norm(fhat.(X) - y,2)^2)
	    push!(g,
	    Axis([
	        Plots.Scatter(X, y, style="mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}"),
	        Plots.Linear(x->fhat([x]), (0,1), style="solid, black, mark=none"),
	        ], title="{\$k=$k, \\text{err} = $(err)\$}", style="xticklabels=\\empty, yticklabels=\\empty", legendPos="outer north east", width="4.5cm", ymin=0, ymax=2))
	end
	g
end

plot(p)
# -

	using Polynomials
	using Distributions
	using QuadGK
	using Random
	using LinearAlgebra

	p = let

		ftrue = x -> x + x^2
		f = x -> ftrue(x) + randn()/2
		Random.seed!(0)
		X = rand(Uniform(0,1), 6)
		y = f.(X)

		k_arr = collect(0:5)
		e_train = Array{Float64}(length(k_arr))
		e_gen = Array{Float64}(length(k_arr))

		for (i,k) in enumerate(k_arr)
			bases = polynomial_bases(1, k)
			B = [b(x) for x in X, b in bases]
    		θ = pinv(B)*y
		    fhat = x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
		    e_train[i] = norm(fhat.(X) - y,2)^2 / length(y)
		    e_gen[i] = QuadGK.quadgk(x->(ftrue(x) - fhat([x]))^2, 0.0, 1.0)[1]
		end

		Axis([
		    Plots.Linear(k_arr[1:5], e_train[1:5], style="solid, pastelBlue, mark=*, mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}", legendentry="training error")
		    Plots.Linear(k_arr, e_gen, style="solid, pastelRed, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="generalization error")
		    ], width="6cm", xlabel="polynomial degree \$k\$", ylabel="error", ymode="log", legendPos="outer north east")
	end

	plot(p)

# +
using LinearAlgebra
using Statistics
import QuadGK: quadgk
p = let

	function regression(X, y, bases)
	    B = [b(x) for x in X, b in bases]
	    θ = pinv(B)*y
	    return x -> sum(θ[i] * bases[i](x) for i in 1 : length(θ))
	end

	polynomial_bases_1d(i, k) = [x->x[i]^p for p in 0:k]
	function polynomial_bases(n, k)
		bases = [polynomial_bases_1d(i, k) for i in 1 : n]
		terms = Function[]
		for ks in Iterators.product([0:k for i in 1:n]...)
			if sum(ks) ≤ k
				push!(terms,
					x->prod(b[j+1](x) for (j,b) in zip(ks,bases)))
			end
		end
		return terms
	end

	xdom = (-5,5)
	f = x -> x/10 + sin(x)/4 + exp(-x^2)

	x_arr = collect(range(-4, stop=4, length=9))
	y_arr = f.(x_arr)

	err_gen(f, fhat, x_arr) = quadgk(x->(f(x) - fhat(x))^2 * 1/(xdom[2] - xdom[1]), xdom...)[1]
	err_train(f, fhat, x_arr) = mean((f.(x_arr) - fhat.(x_arr)).^2)

	G = GroupPlot(3,3,groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, xlabels at=edge bottom, xticklabels at=edge bottom, ylabels at=edge left, yticklabels at=edge left",
	              style="width=4.75cm, ymin=-1, ymax=1.25, enlarge x limits=0, xlabel=\$x\$, ylabel=\$y\$, axis on top=true")

	for k in 0 : 8

	    fhat = regression(x_arr, y_arr, polynomial_bases(1, k))

	    plots = Plots.Plot[]
	    push!(plots, Plots.Linear(f, xdom, style="solid, black, mark=none"))
	    push!(plots, Plots.Scatter(x_arr, y_arr, style="only marks, solid, mark size=1, mark=*, mark options={fill=black, draw=black}"))
	    push!(plots, Plots.Linear(fhat, xdom, style="solid, thick, pastelBlue, mark=none"))
	    push!(plots, Plots.Node("\\footnotesize \$\\epsilon_{train} = $(round(err_train(f, fhat, x_arr), digits=3))\$", 0.4, 0.3, axis="axis description cs", style="right"))
	    push!(plots, Plots.Node("\\footnotesize \$\\epsilon_{gen} = $(round(err_gen(f, fhat, x_arr), digits=3))\$", 0.4, 0.15, axis="axis description cs", style="right"))
	    push!(plots, Plots.Node("\\footnotesize \$k = $k\$", 0.05, 0.87, axis="axis description cs", style="right"))

	    push!(G, Axis(plots))
	end
	G
end

plot(p)

# +
using Random

p = let

	f = x -> sin(x) + randn()*0.1
	Random.seed!(0)
	X = Vector{Float64}[]
	push!(X, [0.1])
	push!(X, [0.2])
	push!(X, [0.35])
	push!(X, [0.4])
	push!(X, [0.6])
	push!(X, [0.8])
	y = (x->f(x[1])).(X)

	fhat = regression(X, y, polynomial_bases(1, 2))

	X2 = Vector{Float64}[]
	push!(X2, [2.6])
	push!(X2, [2.8])
	y2 = (x->f(x[1])).(X2)

	Axis([
	        Plots.Linear(sin, (0,3), style="solid, black, mark=none", legendentry=L"f")
	        Plots.Linear(x->fhat([x]), (0,3), style="solid, pastelBlue, mark=none", legendentry=L"\hat{f}")
	        Plots.Scatter([x[1] for x in X], y, style="only marks, mark=*, mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}", legendentry="training samples")
	        Plots.Scatter([x[1] for x in X2], y2, style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="holdout samples")
	        ], ylabel=L"y", xlabel=L"x", style="xticklabels=\\empty, yticklabels=\\empty, legend cell align=left, legend style={draw=none, at={(0.5,-0.25)}, anchor=north, legend columns=1}", width="9cm")
end

plot(p)

# +
using Random

p = let
	Random.seed!(0)
	xs = range(0, stop=3π, length=10)
	y = sin.(xs)
	X = [[x] for x in xs]

	fit = (X,y) -> regression(X, y, polynomial_bases(1,4))
	metric = (f, X, y) -> begin
	    m = length(X)
	    return sum((f(X[i]) - y[i])^2 for i in m)/m
	end

	sets = [TrainTest(collect(1:9),collect(10:10)),
	         TrainTest(collect(1:8),collect(9:10)),
	         TrainTest(collect(1:5),collect(6:10)),
	        ]

	models = [fit(X[tt.train], y[tt.train]) for tt in sets]
	metrics = [metric(model, X[tt.test], y[tt.test]) for (tt,model) in zip(sets,models)]

	plots = Plots.Plot[]
	push!(plots, Plots.Scatter(xs, y, style="black, mark=*, mark size=1, mark options={draw=black, fill=black}",  legendentry="dataset"))
	push!(plots, Plots.Linear(x->models[1]([x]), (-π/4,3π+π/4), style="solid, mark=none, pastelPurple",  legendentry="10-fold"))
	push!(plots, Plots.Linear(x->models[2]([x]), (-π/4,3π+π/4), style="solid, mark=none, pastelBlue",   legendentry=" 5-fold"))
	push!(plots, Plots.Linear(x->models[3]([x]), (-π/4,3π+π/4), style="solid, mark=none, pastelSeaGreen", legendentry=" 2-fold"))
	Axis(plots, xlabel=L"x", ylabel=L"y", width="10cm", height="4cm", ymin=-5, ymax=10,
	style="legend cell align=left, legend style={draw=none, at={(0.5,-0.15)}, anchor=north, legend columns=4}"
	)
end
plot(p)

# +
using Random

struct TrainTest
	train
	test
end
function k_fold_cross_validation_sets(n, k)
    @assert k ≤ n
    perm = randperm(n)
    sets = TrainTest[]
    for i = 1:k
        validate = perm[i:k:n];
        train = perm[setdiff(1:n, i:k:n)]
        push!(sets, TrainTest(train, validate))
    end
    return sets
end
function train_and_validate(X, y, tt, fit, metric)
    model = fit(X[tt.train], y[tt.train])
    return metric(model, X[tt.test], y[tt.test])
end
function cross_validation_estimate(X, y, sets, fit, metric)
	accuracies = [train_and_validate(X, y, tt, fit, metric) for tt in sets]
    return (mean(accuracies), std(accuracies))
end

p = let
	f = x->sin(2x)*cos(10x)

	Random.seed!(0)
	X = rand(10)
	y = f.(X) + randn(length(X))/10

	sets = k_fold_cross_validation_sets(length(X), 3)
	metric = (f, X, y)->begin
		m = length(X)
		return sum((f(X[i]) - y[i])^2 for i in m)/m
	end

	λ_arr = collect(10 .^ range(-4, stop=2, length=101))
	e_arr = Array{Float64}(undef, length(λ_arr))
	for (i,λ) in enumerate(λ_arr)
	    fit = (X,y)->regression(X, y, radial_bases(r->exp(-5*r^2), X), λ)
	    e_arr[i] = cross_validation_estimate(X, y, sets, fit, metric)[1]
	end

	Axis(Plots.Linear(λ_arr, e_arr, style="solid, black, mark=none"), xmode="log", xlabel=L"\lambda", ylabel="mean cross validated MSE", width="9cm")
end

plot(p)

# +
using Random

function bootstrap_sets(m, b)
	return [TrainTest(rand(1:m, m), 1:m) for i in 1 : b]
end
function bootstrap_estimate(X, y, sets, fit, metric)
	mean(train_and_validate(X, y, tt, fit, metric) for tt in sets)
end
function leave_one_out_bootstrap_estimate(X, y, sets, fit, metric)
	m, b = length(X), length(sets)
	ε = 0.0
	models = [fit(X[tt.train], y[tt.train]) for tt in sets]
	for j in 1 : m
		c = 0
		δ = 0.0
		for i in 1 : b
			if j ∉ sets[i].train
				c += 1
				δ += metric(models[i], [X[j]], [y[j]])
			end
		end
		ε += δ/c
	end
	return ε/m
end
function bootstrap_632_estimate(X, y, sets, fit, metric)
	models = [fit(X[tt.train], y[tt.train]) for tt in sets]
	ϵ_loob = leave_one_out_bootstrap_estimate(X,y,sets,fit,metric)
	ϵ_boot = bootstrap_estimate(X,y,sets,fit,metric)
    return 0.632ϵ_loob + 0.368ϵ_boot
end
function holdout_partition(m, h=div(m,2))
    p = randperm(m)
    train = p[(h+1):m]
    holdout = p[1:h]
    return TrainTest(train, holdout)
end

p = let

	function metric(f, X, y)
	    m = length(X)
	    return sqrt(sum((f(X[i]) - y[i])^2 for i in m)/m)
	end

	Random.seed!(1)
	xs = range(-3, stop=3, length=10)
	f = x -> x^2 + randn()/2
	y = f.(xs)
	X = [[x] for x in xs]

	fit = (X,y) -> regression(X, y, polynomial_bases(1,1))

	mins = Float64[]
	maxs = Float64[]
	medians = Float64[]
	quartile1s = Float64[]
	quartile3s = Float64[]

	bsets() = bootstrap_sets(length(X), 50)
	cvsets() = k_fold_cross_validation_sets(length(X), 5)

	est_methods = [
	    (X,y)->bootstrap_632_estimate(X, y, bsets(), fit, metric),
	    (X,y)->leave_one_out_bootstrap_estimate(X, y, bsets(), fit, metric),
	    (X,y)->bootstrap_estimate(X, y, bsets(), fit, metric),
	    (X,y)->cross_validation_estimate(X, y, cvsets(), fit, metric)[1],
	    (X,y)->train_and_validate(X, y, holdout_partition(length(X), 8), fit, metric),
	]

	for est_method in est_methods
	    vals = [est_method(X,y) for i in 1 : 100]
	    push!(mins, minimum(vals))
	    push!(maxs, maximum(vals))
	    push!(medians, median(vals))
	    push!(quartile1s, quantile(vals, 0.25))
	    push!(quartile3s, quantile(vals, 0.75))
	end

	colors = ["pastelPurple", "pastelBlue", "pastelSeaGreen", "pastelGreen", "pastelRed"]
	plots = Plots.Plot[]
	for i in 1:length(mins)
	    push!(plots, Plots.Command("""\\addplot+[$(colors[i]),
	        boxplot prepared={
	          median=$(medians[i]),
	          upper quartile=$(quartile3s[i]),
	          lower quartile=$(quartile1s[i]),
	          upper whisker=$(maxs[i]),
	          lower whisker=$(mins[i])
	        }, ] coordinates {};"""
	    ))
	end
	Axis(plots, style="ytick={1,2,3,4,5}, yticklabels={0.632 Bootstrap, Leave-One-Out Bootstrap, Bootstrap, 5-Fold Cross Validation, Holdout}", width="8cm", height="4cm", xlabel="root mean squared error")
end
plot(p)

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

# # Sampling Plans
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

using Random
p = let
	L = 6π
	f = x -> (cos(x) - 1)*min(abs(L - x), abs(-x - L))
	X = collect(-L:2π:L)
	Random.seed!(4)
	X2 = X + rand(length(X))*2π .- π
	plots = Plots.Plot[]
	push!(plots, Plots.Linear(f, (-L, L), xbins=151, style="solid, black, mark=none", legendentry=L"f(x)"))
	push!(plots, Plots.Linear(X, f.(X), style="only marks, mark=*, mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}", legendentry="sampling on grid"))
	push!(plots, Plots.Linear(X2, f.(X2), style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="stratified sampling"))
	Axis(plots, width="9cm", xlabel=L"x", ylabel=L"y", style="xtick=\\empty, ytick=\\empty, legend cell align=left, legend style={draw=none, at={(0.5,-0.15)}, anchor=north, legend columns=1}")
end
plot(p)

# +
using Random
using LinearAlgebra

function mutate!(Xs)
    m, n = length(Xs), length(Xs[1])
    j = rand(1:n)
    i = randperm(m)[1:2]
    Xs[i[1]][j], Xs[i[2]][j] = Xs[i[2]][j], Xs[i[1]][j]
    return Xs
end
function uniform_projection_plan(n)
	p = randperm(n)
	[[u,v] for (u,v) in zip(p,1:n)]
end
function pairwise_distances(Xs, p=2)
    n = length(Xs)
    [norm(Xs[i]-Xs[j], p) for i in 1:(n-1) for j in (i+1):n]
end
phiq(Xs, q, p) = sum(d^(-q) for d in pairwise_distances(Xs, p))^(1/q)

p = let

	N = 20
	capacity = 6*3

	Random.seed!(0)
	X = Vector{Vector{Float64}}[]
	push!(X, [Float64[i,i] for i in 1:N])
	for i in 2 : capacity
		push!(X, mutate!(mutate!(deepcopy(X[end]))))
	end
	Φ = [phiq(x, 1, 2) for x in X]

	for i in 1 : 50
		# remove the index with minimum distance to neighbors
		i_min = 1
		d_min = Inf
		p = sortperm(Φ)
		X[:] = X[p]
		Φ[:] = Φ[p]
		for i in 2 : length(X)-1
			d = Φ[i+1] - Φ[i-1]
			@assert d >= 0
			if d < d_min
				d_min = d
				i_min = 1
			end
		end
		X[i_min] = mutate!(X[i_min])
		Φ[i_min] = phiq(X[i_min], 1, 2)
	end

	sort!(X, by=x->phiq(x, 1, 2))

    g = GroupPlot(6,3,groupStyle="horizontal sep=0.5cm, vertical sep=0.75cm")

   	for Y in X
   		push!(g, Axis(
   			Plots.Scatter([x[1] for x in Y], [x[2] for x in Y], style="mark=*, mark size=0.7, mark options={draw=black, fill=black}"),
   			width="4cm", height="4cm", style="ytick=\\empty, xtick=\\empty, title style={yshift=-0.5em}",
   			title = "\$\\Phi_1 = $(round(phiq(Y, 1, 2), digits=1))\$", # $
   			)
		)
   	end

    g
end

plot(p)

# +
using Random

min_dist(a, B, p) = minimum(norm(a-b, p) for b in B)
d_max(A, B, p=2) = maximum(min_dist(a, B, p) for a in A)
function greedy_local_search(X, m, d=d_max)
	S = [X[rand(1:m)]]
	for i in 2 : m
		j = argmin([x ∈ S ? Inf : d(X, push!(copy(S), x))
		           for x in X])
		push!(S, X[j])
	end
	return S
end
function exchange_algorithm(X, m, d=d_max)
	S = X[randperm(length(X))[1:m]]
	δ, done = d(X, S), false
	while !done
		best_pair = (0,0)
		for i in 1 : m
			s = S[i]
			for (j,x) in enumerate(X)
				if !in(x, S)
					S[i] = x
					δ′ = d(X, S)
					if δ′ < δ
						δ = δ′
						best_pair = (i,j)
					end
				end
			end
			S[i] = s
		end
		done = best_pair == (0,0)
		if !done
			i,j = best_pair
			S[i] = X[j]
		end
	end
	return S
end
function multistart_local_search(X, m, alg, k_max, d=d_max)
	sets = [alg(X, m, d) for i in 1 : k_max]
	return sets[argmin([d(X, S) for S in sets])]
end

p = let

	Random.seed!(0)
	X = Vector{Float64}[]
	for i in 1 : 100
		push!(X, rand(2))
	end

	Xs_gls = multistart_local_search(X, 8, greedy_local_search, 10)
	Xs_ea = multistart_local_search(X, 8, exchange_algorithm, 10)

	pX = Plots.Scatter([x[1] for x in X], [x[2] for x in X], style="mark=*, mark size=1, mark options={draw=black, fill=black}")

    g = GroupPlot(2,1,groupStyle="horizontal sep=0.5cm")

    push!(g, Axis(
    	[pX,
    	 Plots.Scatter([x[1] for x in Xs_gls], [x[2] for x in Xs_gls], style="mark=*, mark size=1.5, mark options={draw=pastelRed, fill=pastelRed}")
    	],
    	width="6.25cm", xlabel=L"x_1", ylabel=L"x_2", title="greedy local search",
    	style="ytick=\\empty, xtick=\\empty",
    ))

    push!(g, Axis(
    	[pX,
    	 Plots.Scatter([x[1] for x in Xs_ea], [x[2] for x in Xs_ea], style="mark=*, mark size=1.5, mark options={draw=pastelRed, fill=pastelRed}")
    	],
    	width="6.25cm", xlabel=L"x_1", title="exchange algorithm",
    	style="ytick=\\empty, xtick=\\empty",
    ))


    g
end

plot(p)

# +
using Sobol
import OnlineStats: Mean, Series, value
using Random

p = let

	f = x -> sin(10x)


	Random.seed!(0)
	n_tries_per = 4
	f_true = sin(5)^2 / 5

	arr_n = [round(Int, n) for n in 10 .^ range(0, stop=4, length=101)]
	arr_err_rand = map(arr_n) do n
	    o = Series(Mean())
	    for i in 1 : n_tries_per
	        v = Series(Mean())
	        for j in 1 : n
	            fit!(v, f(rand()))
	        end
	        fit!(o, value(v)[1])
	    end
	    abs(value(o)[1] - f_true) / f_true
	end
	arr_err_sobol = map(arr_n) do n
	    o = Series(Mean())
	    for x in Iterators.take(SobolSeq(1), n)
	        fit!(o, f(x[1]))
	    end
	    abs(value(o)[1] - f_true) / f_true
	end

	Axis([
	    Plots.Linear(arr_n, arr_err_rand, style="mark=none, pastelRed", legendentry="rand"),
	    Plots.Linear(arr_n, arr_err_sobol, style="mark=none, pastelBlue", legendentry="Sobol"),
	    ], xlabel="number of samples", ylabel="relative error", xmode="log", ymode="log", width="9cm",
	       style="trim axis left, trim axis right, legend cell align=left, legend style={draw=none, at={(0.5,-0.4)}, anchor=north, legend columns=2}")
end

plot(p)

# +
using Primes
function halton(i, b)
    result, f = 0.0, 1.0
    while i > 0
        f = f / b;
        result = result + f * mod(i, b)
        i = floor(Int, i / b)
    end
    return result
end
get_filling_set_halton(m; b=2) = [halton(i,b) for i in 1: m]
function get_filling_set_halton(m, n)
    bs = primes(max(ceil(Int, n*(log(n) + log(log(n)))), 6))
    seqs = [get_filling_set_halton(m, b=b) for b in bs[1:n]]
    return [collect(x) for x in zip(seqs...)]
end

p = let

	n = 19
	seqs = [get_filling_set_halton(n, b=b) for b in [19,23]]
	X = [collect(x) for x in zip(seqs...)]

	Axis(
	    Plots.Scatter([x[1] for x in X], [x[2] for x in X],
	                      style="mark=*, mark size=1, mark options={draw=black, fill=black}"),
	    width="5cm", height="5cm",
	style="xtick=\\empty, ytick=\\empty, clip marker paths=true, axis on top=true"
	    )
end

plot(p)

# +
using Primes
using Random
import Printf: @sprintf

function get_filling_set_random_walk(n)
    Xs = [rand()]
    for i in 2 : n
        push!(Xs, mod(Xs[end] + 0.5 + rand(), 1))
    end
    return Xs
end

function get_filling_set_additive_recurrence(n; c=φ-1)
    Xs = [rand()]
    for i in 2 : n
        push!(Xs, mod(Xs[end] + c, 1))
    end
    return Xs
end
function get_filling_set_additive_recurrence(n, k)
    ps = primes(max(ceil(Int, k*(log(k) + log(log(k)))), 6))
    seqs = [get_filling_set_additive_recurrence(n, c=sqrt(p))
            for p in ps[1:k]]
    return [collect(x) for x in zip(seqs...)]
end

p = let

	g = GroupPlot(5,3, groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm")

	ns = [10,100,1000]

	function add_plot(Xs,mark_size)
		push!(g, Axis(
		    Plots.Scatter([x[1] for x in Xs], [x[2] for x in Xs], [i/length(Xs) for i in 1 : length(Xs)],
		                      style="mark=*, mark size=$(mark_size), scatter/use mapped color={draw=mapped color, fill=mapped color}, colormap name=pasteljet"),
		    width="4.5cm", height="4.5cm", enlargelimits=0,
		style="xtick=\\empty, ytick=\\empty, clip marker paths=true, axis on top=true"
		    ))
	end

	Random.seed!(0)
	mark_sizes = [1,0.75,0.5]
	for (j,n) in enumerate(ns)
		mark_size = mark_sizes[j]
	    add_plot([rand(2) for i in 1 : n], mark_size)
	    g.axes[end].ylabel=@sprintf("\\small \$m = %d\$", n)

	    add_plot(get_filling_set_additive_recurrence(n, 2), mark_size)
	    add_plot(get_filling_set_halton(n, 2), mark_size)

	    Xs = Array{Vector{Float64}}(undef, n)
	    for (i,x) in enumerate(Iterators.take(SobolSeq(2), n))
	        Xs[i] = x
	    end
	    add_plot(Xs, mark_size)

	    add_plot(uniform_projection_plan(n)./n, mark_size)
	end

	g.axes[end-4].xlabel="\\small Random"
	g.axes[end-3].xlabel="\\small Additive Recurrence"
	g.axes[end-2].xlabel="\\small Halton"
	g.axes[end-1].xlabel="\\small Sobol"
	g.axes[end].xlabel="\\small Uniform Projection"

	g
end

plot(p)

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

# # Multiobjective
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Vec
using Random

p = let

	Random.seed!(0)
	x_arr = Float64[]
	y_arr = Float64[]
	A = VecE2(1,1)
	for θ in range(180, stop=270, length=90)
	   B = A + polar(0.9 - 0.3rand(), deg2rad(θ))
	   push!(x_arr, B.x)
	   push!(y_arr, B.y)
	end

    Axis([Plots.Command("\\draw[pastelBlue, line width=5mm, line cap=round] (axis cs:0.12,1) [out=270, in=180] to (axis cs:1,0.12);"),
    	  Plots.Linear(x_arr,y_arr, style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}"),
    	  Plots.Node("Ideal",0,0, style="anchor=south west"),
    	  Plots.Node("Pareto Frontier",0.35,0.1, style="anchor=west"),
    	  Plots.Node("Pareto dominated points (suboptimal)",0.45,0.85, style="anchor=west, text width=4cm"),
    	  Plots.Node("criterion space \\\\ (denoted \$\\mathcal{Y}\$) ",0.65,0.6, style="anchor=west, text width=4cm"),
    	], width="10cm", xlabel="alert rate", ylabel="collision rate",
	       xmin=0, xmax=1.1, ymin=0, ymax=1.1,
	       style="xtick=\\empty, ytick=\\empty",
    )
end

plot(p)

# +
using Vec
using Random

function dominates(y, y′)
	lo, hi = extrema(y′ - y)
	return hi > 0 && lo ≥ 0
end
function naive_pareto(xs, ys)
    pareto_xs, pareto_ys = [], []
    for (x,y) in zip(xs,ys)
        if !any(dominates(y′,y) for y′ in ys)
            push!(pareto_xs, x)
            push!(pareto_ys, y)
        end
    end
    return (pareto_xs, pareto_ys)
end

p = let

	Random.seed!(1)
	xs = Vector{Float64}[]
	ys = Vector{Float64}[]
	A = VecE2(1,1)
	for θ in range(180, stop=270, length=90)
	   B = A + polar(0.9 - 0.3rand(), deg2rad(θ))
	   push!(xs, [B.x, B.y])
	   push!(ys, [B.x, B.y])
	end

	pareto, _ = naive_pareto(xs, ys)
	x_arr = [x[1] for x in pareto]
	y_arr = [x[2] for x in pareto]
	p = sortperm(x_arr)
	x_arr[:] = x_arr[p]
	y_arr[:] = y_arr[p]

	filter!(x->x ∉ pareto, xs)


    Axis([Plots.Linear([x[1] for x in xs], [x[2] for x in xs], style="only marks, mark=*, mark size=0.5, mark options={draw=black, fill=black}"),
    	  Plots.Linear(x_arr, y_arr, style="pastelBlue, mark=*, mark size=0.5, mark options={draw=none, fill=pastelBlue}"),
    	], width="9cm", xlabel=L"y_1", ylabel=L"y_2",
	       xmin=0, xmax=1.1, ymin=0, ymax=1.1,
	       style="xtick=\\empty, ytick=\\empty, axis lines=middle",
    )
end

plot(p)
# -

abstract type SelectionMethod end
struct TruncationSelection <: SelectionMethod
	k # top k to keep
end
abstract type MutationMethod end
struct GaussianMutation <: MutationMethod
	σ
end

# +
using Random

function Base.partialsort(t::TruncationSelection, y)
    p = sortperm(y)
    return [p[rand(1:t.k, 2)] for i in y]
end

function rand_population_uniform(M, a, b)
    d = length(a)
    return [a+rand(d).*(b-a) for i in 1:M]
end

abstract type CrossoverMethod end
struct InterpolationCrossover <: CrossoverMethod
    λ
end
crossover(C::InterpolationCrossover, a, b) = (1-C.λ)*a + C.λ*b

function mutate(M::GaussianMutation, child)
    return child + randn(length(child))*M.σ
end

p = let

	f = x -> begin
	    r = 0.5 + 0.5*(2x[1]/(1+x[1]^2))
	    θ = x[2]
	    y1 = 1 - r*cos(θ)
	    y2 = 1 - r*sin(θ)
	    return [y1, y2]
	end
	m_pop = 40
	Random.seed!(0)
	population = rand_population_uniform(m_pop, [-5,-5], [5,5])

	S = TruncationSelection(10)
	C = InterpolationCrossover(0.1)
	M = GaussianMutation(0.1)

	G = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, xticklabels at=edge bottom, yticklabels at=edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$y_1\$, ylabel=\$y_2\$, width=5cm, height=5cm, xmin=0, xmax=2, ymin=0, ymax=2")


	m = length(f(population[1]))
	m_pop = length(population)
	m_subpop = m_pop ÷ m
	for k in 1 : 4
	    fitnesses = f.(population)
	    parents = partialsort(S, [y[1] for y in fitnesses])[1:m_subpop]
	    for i in 2 : m
	       subpop=partialsort(S,[y[i] for y in fitnesses])[1:m_subpop]
	        append!(parents, subpop)
	    end

	    plots = Plots.Plot[]
	    push!(plots, Plots.Linear(y1->1-sqrt(1-(1-y1)^2), (0,1), style="solid, pastelBlue, ultra thick"))
	    push!(plots, Plots.Scatter([y[1] for y in fitnesses],
	                               [y[2] for y in fitnesses], style="clip marker paths, mark=*, mark size=1, mark options={draw=gray, fill=gray}"))
	    push!(G, Axis(plots, title="Gen $k")) #$

	    p = randperm(2m_pop)
	    p_ind=i->parents[mod(p[i]-1,m_pop)+1][(p[i]-1)÷m_pop + 1]
	    parents = [[p_ind(i), p_ind(i+1)] for i in 1 : 2 : 2m_pop]
	    children = [crossover(C,population[p[1]],population[p[2]])
	                    for p in parents]
	    population = [mutate(M, c) for c in children]
	end

	G
end
plot(p)

# +
using Vec
using Random

p = let

	Random.seed!(1)
	ys = Vector{Float64}[]
	A = VecE2(1,1)
	for θ in range(180, stop=270, length=101)
	   B = A + polar(0.9 - 0.3rand(), deg2rad(θ))
	   push!(ys, [B.x, B.y])
	end

	function dominates2(y, y′)
		lo, hi = extrema(y′ - y)
		return hi > 0 && lo ≥ 0
	end
	function get_non_domination_levels(ys)
	    L, m = 0, length(ys)
	    levels = zeros(Int, m)
	    while minimum(levels) == 0
	        L += 1
	        for (i,y) in enumerate(ys)
	            if levels[i] == 0 &&
	               !any((levels[i] == 0 || levels[i] == L) && dominates2(ys[i],y) for i in 1 : m)
	                levels[i] = L
	            end
	        end
	    end
	    return levels
	end

	p = Plots.Plot[]
	levels = get_non_domination_levels(ys)
	for lvl in 1 : maximum(levels)
	    x_arr = [x[1] for x in ys[levels .== lvl]]
	    y_arr = [x[2] for x in ys[levels .== lvl]]
	    color = "pastelBlue!$(round(Int, 100*(maximum(levels) - lvl + 1) / maximum(levels)))" # $
	    push!(p, Plots.Linear(x_arr, y_arr, style="solid, $color, mark=*, mark size=1, mark options={draw=$color, fill=$color}"))
	end
	push!(p, Plots.Node("\\small Level \$1\$",  0.72, 0.20, style="anchor=center"))
	push!(p, Plots.Node("\\small Level \$2\$",  0.99, 0.15, style="anchor=west"))
	push!(p, Plots.Node("\\small Level \$10\$", 0.40, 0.94, style="anchor=west"))
	Axis(p, width="10cm", xlabel=L"y_1", ylabel=L"y_2",
	       xmin=0, xmax=1.2, ymin=0, ymax=1.1,
	       style="xtick=\\empty, ytick=\\empty, axis lines=middle",
	)
end

plot(p)

# +
using Random
using LinearAlgebra

p = let

	rng = MersenneTwister(0)

	function discard_closest_pair!(xs, ys)
		index, min_dist = 0, Inf
		for (i,y) in enumerate(ys)
			for (j, y′) in enumerate(ys[i+1:end])
				dist = norm(y - y′)
				if dist < min_dist
					index, min_dist = rand(rng, [i,j]), dist
				end
			end
		end
		deleteat!(xs, index)
		deleteat!(ys, index)
		return (xs, ys)
	end
	function update_pareto_filter!(filter_xs, filter_ys, xs, ys;
		capacity=length(xs),
		)
	    for (x,y) in zip(xs, ys)
	    	if !any(dominates(y′,y) for y′ in filter_ys)
	            push!(filter_xs, x)
	            push!(filter_ys, y)
	        end
	    end
	    filter_xs, filter_ys = naive_pareto(filter_xs, filter_ys)
	    while length(filter_xs) > capacity
	    	discard_closest_pair!(filter_xs, filter_ys)
	    end
	    return (filter_xs, filter_ys)
	end

	f = x -> begin
	    r = 0.5 + 0.5*(2x[1]/(1+x[1]^2))
	    θ = x[2]
	    y1 = 1 - r*cos(θ)
	    y2 = 1 - r*sin(θ)
	    return [y1, y2]
	end
	m_pop = 40
	Random.seed!(0)
	population = rand_population_uniform(m_pop, [-5,-5], [5,5])

	filter_xs = Vector{Float64}[]
	filter_ys = Vector{Float64}[]

	S = TruncationSelection(10)
	C = InterpolationCrossover(0.1)
	M = GaussianMutation(0.1)

	G = GroupPlot(1,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, xticklabels at=edge bottom, yticklabels at=edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$y_1\$, ylabel=\$y_2\$, width=7cm, height=7cm, xmin=0, xmax=1.5, ymin=0, ymax=1.5")


	K = 10
	m = length(f(population[1]))
	m_pop = length(population)
	m_subpop = m_pop ÷ m
	for k in 1 : K

		fitnesses = f.(population)

		filter_xs, filter_ys = update_pareto_filter!(filter_xs, filter_ys, population, fitnesses)


	    parents = partialsort(S, [y[1] for y in fitnesses])[1:m_subpop]
	    for i in 2 : m
	       subpop=partialsort(S,[y[i] for y in fitnesses])[1:m_subpop]
	        append!(parents, subpop)
	    end

	    if k == K
		    plots = Plots.Plot[]
		    push!(plots, Plots.Linear(y1->1-sqrt(1-(1-y1)^2), (0,1), style="solid, pastelBlue, ultra thick"))
		    push!(plots, Plots.Scatter([y[1] for y in filter_ys],
		                               [y[2] for y in filter_ys], style="clip marker paths, mark=*, mark size=1, mark options={draw=gray, fill=gray}"))
		    push!(G, Axis(plots))
		end

	    p = randperm(2m_pop)
	    p_ind=i->parents[mod(p[i]-1,m_pop)+1][(p[i]-1)÷m_pop + 1]
	    parents = [[p_ind(i), p_ind(i+1)] for i in 1 : 2 : 2m_pop]
	    children = [crossover(C,population[p[1]],population[p[2]])
	                    for p in parents]
	    population = [mutate(M, c) for c in children]
	end

	G
end
plot(p)

# +
using Random

p = let

	rng = MersenneTwister(0)

	function discard_closest_pair!(xs, ys)
		index, min_dist = 0, Inf
		for (i,y) in enumerate(ys)
			for (j, y′) in enumerate(ys[i+1:end])
				dist = norm(y - y′)
				if dist < min_dist
					index, min_dist = rand(rng, [i,j]), dist
				end
			end
		end
		deleteat!(xs, index)
		deleteat!(ys, index)
		return (xs, ys)
	end
	function update_pareto_filter!(filter_xs, filter_ys, xs, ys;
		capacity=length(xs),
		)
	    for (x,y) in zip(xs, ys)
	    	if !any(dominates(y′,y) for y′ in filter_ys)
	            push!(filter_xs, x)
	            push!(filter_ys, y)
	        end
	    end
	    filter_xs, filter_ys = naive_pareto(filter_xs, filter_ys)
	    while length(filter_xs) > capacity
	    	discard_closest_pair!(filter_xs, filter_ys)
	    end
	    return (filter_xs, filter_ys)
	end

	f = x -> begin
	    r = 0.5 + 0.5*(2x[1]/(1+x[1]^2))
	    θ = x[2]
	    y1 = 1 - r*cos(θ)
	    y2 = 1 - r*sin(θ)
	    return [y1, y2]
	end
	m_pop = 40
	Random.seed!(0)
	population = rand_population_uniform(m_pop, [-5,-5], [5,5])

	filter_xs = Vector{Float64}[]
	filter_ys = Vector{Float64}[]

	S = TruncationSelection(10)
	C = InterpolationCrossover(0.1)
	M = GaussianMutation(0.1)

	G = GroupPlot(1,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, xticklabels at=edge bottom, yticklabels at=edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$y_1\$, ylabel=\$y_2\$, width=7cm, height=7cm, xmin=0, xmax=1.5, ymin=0, ymax=1.5")


	dist = 0.15

	K = 10
	m = length(f(population[1]))
	m_pop = length(population)
	m_subpop = m_pop ÷ m
	for k in 1 : K

		fitnesses = f.(population)
		neighbors = zeros(m_pop)
		for (i,y) in enumerate(fitnesses)
			for j in i+1 : m_pop
				y′ = fitnesses[j]
				if norm(y - y′) < dist
					neighbors[i] += 1
					neighbors[j] += 1
				end
			end
		end

		filter_xs, filter_ys = update_pareto_filter!(filter_xs, filter_ys, population, fitnesses)

		for (i,y) in enumerate(fitnesses)
			fitnesses[i] = y*neighbors[i]
		end




	    parents = partialsort(S, [y[1] for y in fitnesses])[1:m_subpop]
	    for i in 2 : m
	        subpop=partialsort(S,[y[i] for y in fitnesses])[1:m_subpop]
	        append!(parents, subpop)
	    end

	    if k == K
		    plots = Plots.Plot[]
		    push!(plots, Plots.Linear(y1->1-sqrt(1-(1-y1)^2), (0,1), style="solid, pastelBlue, ultra thick"))
		    push!(plots, Plots.Scatter([y[1] for y in filter_ys],
		                               [y[2] for y in filter_ys], style="clip marker paths, mark=*, mark size=1, mark options={draw=gray, fill=gray}"))
		    push!(G, Axis(plots))
		end

	    p = randperm(2m_pop)
	    p_ind=i->parents[mod(p[i]-1,m_pop)+1][(p[i]-1)÷m_pop + 1]
	    parents = [[p_ind(i), p_ind(i+1)] for i in 1 : 2 : 2m_pop]
	    children = [crossover(C,population[p[1]],population[p[2]])
	                    for p in parents]
	    population = [mutate(M, c) for c in children]
	end

	G
end
plot(p)

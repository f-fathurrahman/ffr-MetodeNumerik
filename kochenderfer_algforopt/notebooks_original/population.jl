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

# # Population
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

    using Distributions
    using Random
    using LinearAlgebra

    p = let
        function rand_population_uniform(M, a, b)
            d = length(a)
            return [a+rand(d).*(b-a) for i in 1:M]
        end

        using Distributions
        function rand_population_normal(M, μ, Σ)
            D = MvNormal(μ,Σ)
            return [rand(D) for i in 1:M]
        end

        using Distributions
        function rand_population_cauchy(M, μ, σ)
            d = length(μ)
            return [[rand(Cauchy(μ[j],σ[j])) for j in 1:d] for i in 1:M]
        end

        Random.seed!(0)
        m = 1000
        pop1 = rand_population_uniform(m, [-2.0, -2.0], [2.0,2.0])
        pop2 = rand_population_normal(m, [0.0, 0.0], [1.0,1.0])
        pop3 = rand_population_cauchy(m, [0.0, 0.0], [1.0,1.0])
        filter!(x->norm(x) < 10, pop3)

        scatter_style = "clip marker paths, mark=*, mark size=0.75, mark options={draw=pastelBlue, fill=pastelBlue, opacity=0.5}"

        G = GroupPlot(3,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, xticklabels at=edge bottom, yticklabels at=edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$x_1\$, ylabel=\$x_2\$, width=4.75cm, height=4.75cm, xmin=-4, xmax=4, ymin=-4, ymax=4, title style={text height=2ex}")
        push!(G, Axis(
                Plots.Scatter([x[1] for x in pop1], [x[2] for x in pop1], style=scatter_style),
                title = "Uniform",
                ))
        push!(G, Axis(
                Plots.Scatter([x[1] for x in pop2], [x[2] for x in pop2], style=scatter_style),
                title = "Normal",
                ))
        push!(G, Axis(
                Plots.Scatter([x[1] for x in pop3], [x[2] for x in pop3], style=scatter_style),
                title = "Cauchy",
                ))
        G
    end

    plot(p)

    using Random
    using Distributions

    function rand_population_uniform(M, a, b)
        d = length(a)
        return [a+rand(d).*(b-a) for i in 1:M]
    end

    abstract type SelectionMethod end
    struct TruncationSelection <: SelectionMethod
        k # top k to keep
    end
    function select(t::TruncationSelection, y)
        p = sortperm(y)
        return [p[rand(1:t.k, 2)] for i in y]
    end

    abstract type CrossoverMethod end
    struct SinglePointCrossover <: CrossoverMethod end
    function crossover(::SinglePointCrossover, a, b)
        i = rand(1:length(a))
        return vcat(a[1:i], b[i+1:end])
    end

    abstract type MutationMethod end
    struct GaussianMutation <: MutationMethod
        σ
    end
    mutate(U::GaussianMutation, child) = child + randn(length(child))*U.σ

    p = let
        K = 4
        G = GroupPlot(K,1,groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, ylabels at=edge left, xlabels at=edge bottom",
                          style="width=5cm, height=5cm, xlabel=\$x_1\$, ylabel=\$x_2\$, xtick=\\empty, ytick=\\empty, contour/labels=false, view={0}{90}")

        function michalewicz(x; m=10)
            return -sum(sin(v)*sin(i*v^2/π)^(2m) for
                       (i,v) in enumerate(x))
        end
        f = x -> michalewicz(x)
        xdomain = (0,4)
        ydomain = (0,4)

        S = TruncationSelection(10)
        C = SinglePointCrossover()
        M = GaussianMutation(0.1)

        m = 40
        Random.seed!(0)
        population = rand_population_uniform(m, [0.0, 0.0], [4.0,4.0])

        for i in 1 : K

            push!(G, Axis([
                Plots.Image((x,y)->f([x,y]), xdomain, ydomain, xbins=600, ybins=600, colormap = pasteljet, colorbar = false),
                Plots.Scatter([P[1] for P in population], [P[2] for P in population], style="mark=*, mark size=1, mark options={draw=black, fill=black}"),
                ]))

            parents = select(S, f.(population))
            children = [crossover(C,population[p[1]],population[p[2]])
                        for p in parents]
            population = [mutate(M, c) for c in children]
        end

        G
    end

    plot(p)

    using Random
    using StatsBase

    p = let
        N = 4
        g = GroupPlot(N,2,groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, ylabels at=edge left, xlabels at=edge bottom",
                          style="width=5cm, height=5cm, xlabel=\$x_1\$, ylabel=\$x_2\$, xtick=\\empty, ytick=\\empty, contour/labels=false, view={0}{90}")

        f = x -> -20*exp(-0.2*sqrt(0.5)*norm(x)) - exp(0.5*(cos(2π*x[1]) + cos(2π*x[2]))) + ℯ + 20
        f2 = (x,y) -> f([x,y])
        xdomain = (-5,5)
        ydomain = (-5,5)
        P = 0.5
        w = 0.2

        function tick_differential_evolution!(f, population, P, w)
            n, m = length(population[1]), length(population)
            for (k,x) in enumerate(population)
                weights = Weights([j!=k for j in 1 : m])
                a, b, c = sample(population, weights, 3, replace=false)
                z = a + w*(b-c)
                j = rand(1:n)
                x′ = [i == j || rand() < P ? z[i] : x[i] for i in 1:n]
                if f(x′) < f(x)
                    x[:] = x′
                end
            end
            return population
        end

        Random.seed!(0)
        population = [(rand(2) .- 0.5).*9 for i in 1 : 20]

        for i in 1 : 2N

            ax = Axis([
                Plots.Image(f2, xdomain, ydomain, xbins=600, ybins=600, colormap = pasteljet, colorbar = false),
                Plots.Scatter([P[1] for P in population], [P[2] for P in population], style="mark=*, mark size=1, mark options={draw=black, fill=black}"),
                ])
            push!(g, ax)

            tick_differential_evolution!(f, population, P, w)
        end

        g
    end

    plot(p)

    using Random
    using Vec
    struct Particle
        x
        v
        x_best
    end
    function initialize_particles_uniform(f, N, a_low, a_high)
        m = length(a_low)
        Δ = a_high - a_low
        spawn = () -> begin
            x = a_low + rand(m).*Δ
            v = (rand(m) .- 0.5) .* 2Δ
            Particle(x, v, x)
        end
        [spawn() for i in 1 : N]
    end

    p = let

        g = GroupPlot(4,2,groupStyle="xlabels at=edge bottom, ylabels at =edge left, horizontal sep=0.5cm, vertical sep=0.5cm")

        f = x -> -exp(-(x[1]*x[2] - 1.5)^2 -(x[2]-1.5)^2)

        xdomain = (0, 3)
        ydomain = (0, 3)

        p_contour = Plots.Contour(f, xdomain, ydomain, style="width=\\columnwidth")

        function add_particles!(particles)

            plots = Plots.Plot[]
            push!(plots, p_contour)
            push!(plots, Plots.Scatter([P.x[1] for P in particles], [P.x[2] for P in particles], style="mark=*, mark size=1, mark options={draw=black, fill=black}"))
            for P in particles
                A = VecE2(P.x[1], P.x[2])
                B = A + VecE2(P.v[1], P.v[2])/8
                push!(plots, Plots.Linear([A.x, B.x], [A.y, B.y], style="solid, black, mark=none"))
            end

            ax = Axis(plots, ymin=0, xmin=0, ymax=3, xmax=3, width="5cm", height="5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
            push!(g, ax)
        end

        Random.seed!(0)
        N = 40
        population = initialize_particles_uniform(f, N, [xdomain[1],ydomain[1]], [xdomain[2],ydomain[2]])
        add_particles!(population)

        w, c1, c2, m = 0.1, 0.25, 2.0, 2
        x_best, y_best = copy(population[1].x_best), Inf
        for P in population
            y = f(P.x)
            if y < y_best; x_best[:], y_best = P.x, y; end
        end
        for iter in 1 : g.dimensions[1]*g.dimensions[2]-1

            for P in population
                r1, r2 = rand(m), rand(m)
                P.x[:] = P.x + P.v
                P.v[:] = w*P.v + c1*r1.*(P.x_best - P.x) + c2*r2.*(x_best - P.x)
                y = f(P.x)
                if y < y_best; x_best[:], y_best = P.x, y; end
                if y < f(P.x_best); P.x_best[:] = P.x; end
            end
            add_particles!(population)
        end

        g
    end

    plot(p)

    using Distributions
    using Random

    p = let
        G = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, horizontal sep=0.5cm, vertical sep=0.5cm")

        function branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π))
            return a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
        end
        f = x -> branin(x)

        xdomain = (-5, 10)
        ydomain = ( 0, 15)

        p_contour = Plots.Contour(f, xdomain, ydomain, levels=[1,2,3,5,10,20,50,100], xbins=101, ybins=101, style="width=\\columnwidth")

        function add_axis!(population)
            plots = Plots.Plot[]
            push!(plots, p_contour)
            push!(plots, Plots.Scatter([x[1] for x in population], [x[2] for x in population], style="mark=*, mark size=1, mark options={draw=black, fill=black}"))
            ax = Axis(plots, ymin=ydomain[1], xmin=xdomain[1], ymax=ydomain[2], xmax=xdomain[2], width="5cm", height="5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
            push!(G, ax)
        end

        Random.seed!(0)
        population = [[rand(Uniform(xdomain[1], xdomain[2])),
                       rand(Uniform(ydomain[1], ydomain[2]))] for i in 1 : 40]
        add_axis!(population)


        α = 0.5
        β = 1.0
        I = r->exp(-0.1r^2)

        N = MvNormal(diagm(0=>ones(length(population[1]))))
        for K in 1 : 3
            for a in population, b in population
                if f(b) < f(a)
                    r = norm(b-a)
                    a[:] += β*I(r)*(b-a) + α*rand(N)
                end
            end
            add_axis!(population)
        end

        G
    end

    plot(p)

    using Distributions
    using Random

    mutable struct Nest
        x # position
        y # value, f(x)
    end

    p = let
        G = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, horizontal sep=0.5cm, vertical sep=0.5cm")

        function branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π))
            return a*(x[2]-b*x[1]^2+c*x[1]-r)^2 + s*(1-t)*cos(x[1]) + s
        end
        f = x -> branin(x)

        xdomain = (-5, 10)
        ydomain = ( 0, 15)

        p_contour = Plots.Contour(f, xdomain, ydomain, levels=[1,2,3,5,10,20,50,100], xbins=101, ybins=101, style="width=\\columnwidth")

        function add_axis!(nests)
            plots = Plots.Plot[]
            push!(plots, p_contour)
            push!(plots, Plots.Scatter([P.x[1] for P in nests], [P.x[2] for P in nests], style="mark=*, mark size=1, mark options={draw=black, fill=black}"))
            ax = Axis(plots, ymin=ydomain[1], xmin=xdomain[1], ymax=ydomain[2], xmax=xdomain[2], width="5cm", height="5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
            push!(G, ax)
        end

        Random.seed!(0)
        N = 40
        population = Array{Nest}(undef, N)
        for i in 1 : N
            x = [xdomain[1] + rand()*(xdomain[2] - xdomain[1]),
                 ydomain[1] + rand()*(ydomain[2] - ydomain[1])]
            population[i] = Nest(x, f(x))
        end
        add_axis!(population)

        a = 30
        s = 0.5
        C = Distributions.Cauchy(0,1.0)

        N, m = length(population), length(population[1].x)
        for iter in 1 : 3
            # get a cuckoo randomly and generate a new nest by Levy flight
            # replace another random cuckoo with that new position if it is better
            i, j = rand(1:N), rand(1:N)
            x = population[j].x + s*[rand(C) for k in 1 : m]
            y = f(x)
            if y < population[i].y
                population[i].x[:] = x
                population[i].y = y
            end

            p = sortperm(population, by=nest->nest.y, rev=true)
            for i in 1 : a
                j = rand(1:N-a)+a
                population[p[i]].x = population[p[j]].x +
                                     s*[rand(C) for k in 1 : m]
                population[p[i]].y = f(population[p[i]].x)
            end
            add_axis!(population)
        end

        G
    end

    plot(p)

        p = let

            G = GroupPlot(2,1,groupStyle="horizontal sep=1cm")

            f = x -> -exp(-x^2) + -2exp(-(x-3)^2)
            x = [-1.8, -1, -0.4, 0.6, 1.0, 1.2]
            x2 = x .* 0.3

            p = Plots.Plot[]
            push!(p, Plots.Linear(f, (-2,5), style="solid, black, mark=none"))
            push!(p, Plots.Scatter(x2, f.(x2), style="black, mark=*, mark size=1.5, mark options={draw=black, fill=black}"))
            push!(G, Axis(p, xlabel=L"x", ylabel=L"y", style="enlarge x limits=0", width="6cm", title="Lamarckian"))

            p = Plots.Plot[]
            push!(p, Plots.Linear(f, (-2,5), style="solid, black, mark=none"))
            for i in 1 : length(x)
               push!(p, Plots.Linear([x[i], x[i]], [f(x[i]), f(x2[i])], style="solid, pastelBlue, mark=*, mark size=1.45, mark options={draw=pastelBlue, fill=white}"))
            end
            push!(p, Plots.Scatter(x, f.(x), style="black, mark=*, mark size=1.5, mark options={draw=black, fill=black}"))

            push!(G, Axis(p, xlabel=L"x", style="enlarge x limits=0", width="6cm", title="Baldwinian"))
        end
        plot(p)

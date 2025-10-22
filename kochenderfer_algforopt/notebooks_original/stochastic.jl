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

# # Stochastic
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Random
using LinearAlgebra

p = let

	xdomain = (-1.6, 1.6)
	ydomain = (-1.6, 1.6)

	A = [1 0; 0 -1]

	f = x -> dot(x, A*x)
	∇f = x -> (A + A')*x

	Random.seed!(2)
	x0 = Float64[1.5, 0.0]
	pts_SD = Vector{Float64}[x0]
	pts_SGD = Vector{Float64}[x0]
	α = 0.2
	σ = 0.1
	for i in 1 : 8
	    x = pts_SD[end]
	    g = ∇f(x)
	    push!(pts_SD, x - α*g)

	    x = pts_SGD[end]
	    g = ∇f(x)
	    push!(pts_SGD, x - α*g + σ*randn(2))
	end

	plots = Plots.Plot[]
	# levels=[-1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2]
	push!(plots, Plots.Contour(f, xdomain, ydomain, style="width=\\columnwidth, forget plot", xbins=101, ybins=101))
	push!(plots, Plots.Linear3([p[1] for p in pts_SGD],
	                           [p[2] for p in pts_SGD],
	                           [f(p) for p in pts_SGD], style="pastelRed, ultra thick, solid, mark=none", legendentry="stochastic gradient descent")),
	push!(plots, Plots.Linear3([p[1] for p in pts_SD],
	                           [p[2] for p in pts_SD],
	                           [f(p) for p in pts_SD], style="pastelBlue, ultra thick, solid, mark=none", legendentry="steepest descent"))


	Axis(plots, width="9cm", height="9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}, legend cell align=left, legend style={draw=none, at={(0.5,-0.15)}, anchor=north, legend columns=1},")
end

plot(p)

# +
using Random
using LinearAlgebra

struct Pattern
    center::Vector{Float64}
    fringe::Vector{Vector{Float64}}
    sidelen::Float64
end

p = let

	function rand_positive_spanning_set(α, n)
	    δ = round(Int, 1/sqrt(α))
	    L = Matrix(Diagonal(δ*rand([1,-1], n)))
	    for i in 1 : n-1
	    	for j in i+1:n
	    		L[i,j] = rand(-δ+1:δ-1)
	    	end
	    end
	    D = L[randperm(n),:]
	    D = L[:,randperm(n)]
	    D = hcat(D, -sum(D,dims=2))
	    return [D[:,i] for i in 1 : n+1]
	end

	g = GroupPlot(4,2,groupStyle="xlabels at=edge bottom, ylabels at=edge left, horizontal sep=0.25cm, vertical sep=0.25cm")

	f = x -> -exp(-(x[1]*x[2] - 1.5)^2 -(x[2]-1.5)^2)

	xdomain = (0, 3)
	ydomain = (0, 3)

	p_contour = Plots.Contour(f, xdomain, ydomain, xbins=151, ybins=151)

	function add_pattern!(patterns)

	    plots = Plots.Plot[]
	    push!(plots, p_contour)

	    for (i,P) in enumerate(patterns)
	        for c in P.fringe
	            push!(plots, Plots.Linear3([P.center[1], c[1]],
	                                       [P.center[2], c[2]],
	                                       [f(P.center), f(c)],
	                style="solid, thick, black, mark=*, mark size=1, opacity=$(0.3^(length(patterns) - i)), mark options={solid, draw=black}"))
	        end
	    end

	    ax = Axis(plots, ymin=0, xmin=0, ymax=3, xmax=3, width="4.7cm", height="4.7cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
	    push!(g, ax)
	end

	x = [1.75,1.25]
	α, y, n = 1, f(x), length(x)
	patterns = Pattern[]

	Random.seed!(4)
	for iter in 1 : g.dimensions[1]*g.dimensions[2]
	    improved = false
	    D = rand_positive_spanning_set(α, n)
	    push!(patterns, Pattern(x, [x + α*d for d in D], sqrt(α)))
	    add_pattern!(patterns)
	    for (i,d) in enumerate(D)
	        x′ = x + α*d
	        y′ = f(x′)
	        if y′ < y
	            x, y, improved = x′, y′, true

	            x′ = x + 3α*d
	            y′ = f(x′)
	            if y′ < y
	                x, y = x′, y′
	            end
	            break
	        end
	    end
	    α = improved ? min(4α, 1) : α/4
	end

	g
end

plot(p)
# -

            using Distributions
            p = let
                Axis([Plots.Linear(x->pdf(Exponential(0.5),x), (0,5), style="thick, pastelSeaGreen", legendentry=L"\lambda = 1/2"),
                      Plots.Linear(x->pdf(Exponential(1.0),x), (0,5), style="thick, pastelBlue", legendentry=L"\lambda = 1"),
                      Plots.Linear(x->pdf(Exponential(2.0),x), (0,5), style="thick, pastelPurple", legendentry=L"\lambda = 2"),
                    ], xlabel=L"x", ylabel=L"\lambda e^{-\lambda x}", ymin=0, width="9cm", style="enlarge x limits=0, legend style={draw=none, at={(0.5,-0.25)}, anchor=north, legend columns=1}")
            end
            plot(p)

# +
p = let
t_arr = collect(10 .^ range(0.0, stop=4.0, length=101))
T₀ = 10.0
T1 = x->T₀*log(2)/log(x+1) # logarithmic annealing
T2 = (x,γ)->T₀*γ^(x-1) # exp annealing
T3 = x->T₀/x # fast annealing
Axis([
Plots.Linear(t_arr, T1.(t_arr), style="solid, pastelRed, mark=none", legendentry="logarithmic"),
Plots.Linear(t_arr, T2.(t_arr, 0.25), style="solid, pastelBlue, mark=none", legendentry="exponential, \$\\gamma=1/4\$"),
Plots.Linear(t_arr, T2.(t_arr, 0.50), style="solid, pastelBlue!70, mark=none", legendentry="exponential, \$\\gamma=1/2\$"),
Plots.Linear(t_arr, T2.(t_arr, 0.75), style="solid, pastelBlue!40, mark=none", legendentry="exponential, \$\\gamma=3/4\$"),
Plots.Linear(t_arr, T3.(t_arr), style="solid, pastelGreen, mark=none", legendentry="fast"),
], xmode="log", style="axis on top, enlarge x limits=0, legend cell align=left, legend style={draw=none, at={(0.5,-0.35)}, anchor=north, legend columns=1,},",
xlabel="iteration", ylabel="temperature", width="1.1*9cm",
)
end

plot(p)

# +
using Distributions
using Random

struct OnlineMeanAndVariance
	count::Int
	mean::Float64
	M2::Float64
end

function update(omav::OnlineMeanAndVariance, v::Float64)
    count, mean, M2 = omav.count, omav.mean, omav.M2
    count += 1
    delta = v - mean
    mean += delta / count
    delta2 = v - mean
    M2 += delta * delta2

    return OnlineMeanAndVariance(count, mean, M2)
end

get_μ(omav::OnlineMeanAndVariance) = omav.mean
get_ν(omav::OnlineMeanAndVariance) = omav.M2/omav.count

p = let

	function ackley(x, a=20, b=0.2, c=2π)
	    d = length(x)
		return -a*exp(-b*sqrt(sum(x.^2)/d)) -
		          exp(sum(cos.(c*xi) for xi in x)/d) + a + ℯ
	end

	function simulated_annealing!(f, x, P, t, y_arr)
	    y = f(x)
	    x_best, y_best = copy(x), y
	    y_arr[1] = y_best
	    for k in 1 : length(y_arr)-1
	        x′ = x + rand(P)
	        y′ = f(x′)
	        Δy = y′ - y
	        if Δy ≤ 0 || rand() < exp(-Δy/t(k))
	            x, y = x′, y′
	        end
	        if y′ < y_best
	            x_best, y_best = x′, y′
	        end
	        y_arr[1+k] = y_best
	    end
	    return x_best
	end

	Ps = [MvNormal([1.0,1.0]), MvNormal([5.0,5.0]), MvNormal([25.0,25.0])]
	ts = [k->1/k, k->10/k, k->25/k]

	x0 = [15.0,15.0]
	N = 100

	ntrials = 10
	μν_arr = Vector{OnlineMeanAndVariance}(undef, N+1)
	y_arr = Array{Float64}(undef, N+1)
	k_arr = collect(1:N+1)

	Ps = [MvNormal([1.0,1.0]), MvNormal([5.0,5.0]), MvNormal([25.0,25.0])]
	P_titles = [1, 5, 25]
	ts = [k->1/k, k->10/k, k->25/k]
	t0s = [1,10,25]

	G = GroupPlot(3,3,groupStyle="horizontal sep=0.25cm, vertical sep=1cm, xlabels at=edge bottom, xticklabels at=edge bottom, ylabels at=edge left, yticklabels at=edge left", style="width=4.25cm, xlabel=iteration, ylabel=\$y\$, ymin=-5, ymax=30, xmin=1, xmax=100")

	Random.seed!(0)
	for (i,(P,σ)) in enumerate(zip(Ps,P_titles))
	    for (j,(t,t0)) in enumerate(zip(ts, t0s))
	        for s in 1 : N+1
	            μν_arr[s] = OnlineMeanAndVariance(0,0.0,0.0)
	        end
	        for k in 1 : ntrials
	            simulated_annealing!(ackley, copy(x0), P, t, y_arr)
	            for (s,v) in enumerate(y_arr)
	                μν_arr[s] = update(μν_arr[s], v)
	            end
	        end

	        mean_arr = [get_μ(omav) for omav in μν_arr]
	        q05_arr = [quantile(Normal(get_μ(omav), sqrt(get_ν(omav))+0.01), 0.05) for omav in μν_arr]
	        q25_arr = [quantile(Normal(get_μ(omav), sqrt(get_ν(omav))+0.01), 0.25) for omav in μν_arr]
	        q75_arr = [quantile(Normal(get_μ(omav), sqrt(get_ν(omav))+0.01), 0.75) for omav in μν_arr]
	        q95_arr = [quantile(Normal(get_μ(omav), sqrt(get_ν(omav))+0.01), 0.95) for omav in μν_arr]

	        plots = Plots.Plot[]
	        push!(plots, Plots.Linear(k_arr, q05_arr, style="name path=A, draw=none, mark=none"))
	        push!(plots, Plots.Linear(k_arr, q95_arr, style="name path=B, draw=none, mark=none"))
	        push!(plots, Plots.Command(L"\addplot [pastelBlue!50, opacity=0.5] fill between[of=A and B];"))
	        push!(plots, Plots.Linear(k_arr, q25_arr, style="name path=A, draw=none, mark=none"))
	        push!(plots, Plots.Linear(k_arr, q75_arr, style="name path=B, draw=none, mark=none"))
	        push!(plots, Plots.Command(L"\addplot [pastelBlue!50, opacity=0.5] fill between[of=A and B];"))
	        push!(plots, Plots.Linear(k_arr, mean_arr, style="solid, black, mark=none"))
	        push!(G, Axis(plots, title="\$\\sigma = $σ, t^{(1)} = $(t0)\$"))
	    end
	end

	G
end

plot(p)
# -

    p = let

        c = 2
        f = x -> (1 + c*(x-0.6)/0.4)*(x > 0.6) +
                 (1 + c*(0.4-x)/0.4)^(-1)*(x < 0.4) +
                 (0.4 ≤ x ≤ 0.6)

        Axis(Plots.Linear(f, (0,1), style="solid, pastelBlue, mark=none"),
             xlabel="Acceptance Ratio \$a_i / n_s\$",
             ylabel="Step Multiplication Factor",
             ymin=0,
             xmin=0, xmax=1,
             width="1.1*9cm",
             style = "xtick={0,0.4,0.6,1}, ytick={"*string(f(0))*",1,"*string(f(1))*"}, yticklabels={\$\\frac{1}{1+c}\$,1,\$c\$}")
    end

    plot(p)

	using Random
    using Distributions
    using Discretizers

    struct SimAnnealParticle
        x::Float64
        x_best::Float64
    end

    p = let

        f = x -> (6 - 5exp(-x^2) + sin(10x))/3

        function update_particle(X::SimAnnealParticle, P::Normal, t::Float64)
            x, y = X.x, f(X.x)
            x_best, y_best = X.x_best, f(X.x_best)
            x′ = x + rand(P)
            y′ = f(x′)
            Δy = y′ - y
            if Δy ≤ 0 || rand() < exp(-Δy/t)
                x, y = x′, y′
            end
            if y′ < y_best
                x_best, y_best = x′, y′
            end
            return SimAnnealParticle(x, x_best)
        end
        function construct_histogram_linear_data(
            data::Vector{Q},
            binedges::Vector{R},
            ) where {Q<:Real, R<:Real}

            n = length(binedges)
            disc = LinearDiscretizer(binedges)
            counts = get_discretization_counts(disc, data)

            arr_x = convert(Vector{Float64}, binedges)
            arr_y = convert(Vector{Float64}, counts)
            arr_y ./= sum(counts)
            arr_y ./= binwidths(disc)
            push!(arr_y, arr_y[end])

            Plots.Linear(hcat(arr_x, arr_y)', style="ybar interval, fill=pastelBlue!50, draw=none", mark="none")
        end

        G = GroupPlot(4,2, groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, xlabels at=edge bottom, xticklabels at=edge bottom, ylabels at=edge left, yticklabels at=edge left", style="width=5cm, xlabel=\$x\$, ylabel=\$y\$, ymin=0, ymax=2.5, xmin=-2, xmax=2, xtick=\\empty, ytick=\\empty, axis on top")

        nbins = 100
        nparticles = 10000
        nsteps = 8

        P = Normal(0.25,1.0)
        t = 1.0
        t_decay = 0.5

        particles = [SimAnnealParticle(-1.5, -1.5) for i in 1 : nparticles]

        Random.seed!(0)
        for k in 1 : nsteps

            plots = Plots.Plot[]
            push!(plots, construct_histogram_linear_data([P.x for P in particles if abs(P.x) < 2], collect(range(-2, stop=2, length=nbins+1))))
            push!(plots, Plots.Linear(f, (-2, 2), style="solid, thick, black, mark=none"))
            push!(G, Axis(plots))

            for (i, X) in enumerate(particles)
                particles[i] = update_particle(X, P, t)
            end
            t *= t_decay
        end

        G
    end

    plot(p)

  using Distributions
  using Random

  p = let

      g = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$x_1\$, ylabel=\$x_2\$")

      branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π)) = a*(x[2] - b*x[1]^2 + c*x[1] - r)^2 + s*(1-t)*cos(x[1]) + s
      f = x -> branin(x)
      xdomain = (-5, 10)
      ydomain = ( 0, 15)

      p_image = Plots.Image((x,y)->f([x,y]), xdomain, ydomain, colormap=pasteljet, colorbar=false, xbins=600, ybins=600)

      function add_frame!(samples, P, p, elite_size)
      	  elites = p[1:elite_size]
      	  plebians = p[elite_size+1 : end]

          plots = Plots.Plot[]
          push!(plots, p_image)
          push!(plots, Plots.Contour(x->pdf(P, x), xdomain, ydomain, contour_style="draw color=white", style="white", xbins=151, ybins=151))
          push!(plots, Plots.Scatter(samples[1,plebians], samples[2,plebians], style="clip marker paths, mark=*, mark size=0.75, mark options={draw=white, fill=white}"))
          push!(plots, Plots.Scatter(samples[1,elites], samples[2,elites], style="clip marker paths, mark=*, mark size=0.75, mark options={draw=pastelRed, fill=pastelRed}"))

          ax = Axis(plots, ymin=ydomain[1], xmin=xdomain[1], ymax=ydomain[2], xmax=xdomain[2], width="5cm", height="5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
          push!(g, ax)
      end

      Random.seed!(0)
      N = 40
      elite_size = 10
      P = MvNormal([mean(xdomain), mean(ydomain)], Matrix{Float64}(10I, 2, 2))
      for iter in 1 : g.dimensions[1]*g.dimensions[2]
          samples = rand(P, N)
          p = sortperm([f(samples[:,j]) for j in 1:N])
          add_frame!(samples, P, p, elite_size)
          P = fit(MvNormal, samples[:,p[1:elite_size]])
      end

      g
  end

  plot(p)

    using Distributions
    using Random

    p = let
        xdom = (-6,6)
        f = x -> (1.1 - exp(-(x-3)^2) - exp(-(x+3)^2))/4

        Random.seed!(3)
        GMM = MixtureModel([Normal(-3.0,1.0), Normal(3.0,1.0)])
        N = fit_mle(Normal, rand(GMM, 1000))
        samples = rand(GMM, 10)

        plots = Plots.Plot[]
        push!(plots, Plots.Linear(f, xdom, style="solid, thick, black, mark=none", legendentry=L"f"))
        push!(plots, Plots.Linear(x->pdf(GMM, x), xdom, style="solid, pastelBlue, mark=none", legendentry="mixture model fit"))
        push!(plots, Plots.Linear(x->pdf(N, x), xdom, style="solid, pastelRed, mark=none", legendentry="normal fit"))
        push!(plots, Plots.Scatter(samples, f.(samples), style="only marks, black, mark=*, mark size=1", legendentry="samples"))
        Axis(plots, width="6cm", xmin=xdom[1], xmax=xdom[2], ymin=0, xlabel=L"x", style="legend pos = outer north east")
    end

    plot(p)

    using Distributions
    using Random

    p = let
        wheeler(x, a=1.5) = - exp(-(x[1]*x[2] - a)^2 -(x[2]-a)^2)
        f = x -> wheeler(x)
        xdomain = (0,3)
        ydomain = (0,3)

        G = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at =edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$x_1\$, ylabel=\$x_2\$")

        function add_axis!(samples, P)
            plots = Plots.Plot[]
            push!(plots, Plots.Image((x,y)->f([x,y]), xdomain, ydomain, colormap=pasteljet, colorbar=false, xbins=600, ybins=600))
            push!(plots, Plots.Contour(x->pdf(P, x), xdomain, ydomain, contour_style="draw color=white", style="white", xbins=151, ybins=151))
            push!(plots, Plots.Scatter([x[1] for x in samples], [x[2] for x in samples], style="mark=*, mark size=0.5, mark options={draw=white, fill=white}"))
            ax = Axis(plots, ymin=ydomain[1], xmin=xdomain[1], ymax=ydomain[2], xmax=xdomain[2], width="5cm", height="5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
            push!(G, ax)
        end

        ∇logP_μ(x, μ, Σ⁻¹) = Σ⁻¹*(x - μ)
        ∇logP_Σ(x, μ, Σ⁻¹) = (Σ⁻¹*(x-μ)*(x-μ)'*Σ⁻¹ - Σ⁻¹).*0.5
        function ∇logP_A(x, μ, A, Σ⁻¹)
            ∇Σ = ∇logP_Σ(x, μ, Σ⁻¹)
            return A*(∇Σ + ∇Σ')
        end

        Random.seed!(0)
        A = Matrix(0.5I, 2, 2)
        Σ = A'*A
        σ = 1.0
        μ = [2.0, 2.0] # initial point
        α = 1.0
        M = 40 # number of samples

        for k in 1 : prod(G.dimensions)
            Σ⁻¹ = inv(Σ + Matrix(1.0I, 2, 2)./5)
            P = MvNormal(μ, Σ)
            x_arr = [rand(P) for i in 1 : M]
            f_arr = f.(x_arr)
            add_axis!(x_arr, P)
            μ -= α * sum(v*∇logP_μ(x, μ, Σ⁻¹) for (v,x) in zip(f_arr, x_arr))/M
            A -= α * sum(v*∇logP_A(x, μ, A, Σ⁻¹) for (v,x) in zip(f_arr, x_arr))/M
            Σ = A'*A
        end

        G
    end

    plot(p)

# +
using Distributions
using Random

p = let
	wheeler(x, a=1.5) = - exp(-(x[1]*x[2] - a)^2 -(x[2]-a)^2)
	f = x -> wheeler(x)
	xdomain = (0,3)
	ydomain = (0,3)

	Random.seed!(0)

	μ = [1.75,1.75]
	σ = 0.7
	Σ = [1.0 0.0; 0.0 1.0]

	m = 4 + floor(Int, 3*log(length(μ)))
	m_elite = div(m,2)
	n = length(μ)

	P = MvNormal(μ, σ^2*Σ)
	xs = [rand(P) for i in 1 : m]
	ys = f.(xs)
	is = sortperm(ys) # best to worst

	ws = log((m+1)/2) .- log.(1:m)
	ws[1:m_elite] ./= sum(ws[1:m_elite])

	μ_eff = 1 / sum(ws[1:m_elite].^2)
	cσ = (μ_eff + 2)/(n + μ_eff + 5)
	dσ = 1 + 2max(0, sqrt((μ_eff-1)/(n+1))-1) + cσ
	cΣ = (4 + μ_eff/n)/(n + 4 + 2μ_eff/n)
	c1 = 2/((n+1.3)^2 + μ_eff)
	cμ = min(1-c1, 2*(μ_eff-2+1/μ_eff)/((n+2)^2 + μ_eff))
	ws[m_elite+1:end] .*= -(1 + c1/cμ)/sum(ws[m_elite+1:end])

	δs = [(x - μ)/σ for x in xs]
	δw = sum(ws[i]*δs[is[i]] for i in 1 : m_elite)
	μ_cma = μ + σ*δw

	μ_ce = mean(xs[is[i]] for i in 1:m_elite)

	plots = Plots.Plot[]
	push!(plots, Plots.Image((x,y)->f([x,y]), xdomain, ydomain, colormap=pasteljet, colorbar=false, xbins=600, ybins=600))
	push!(plots, Plots.Contour(x->pdf(P, x), xdomain, ydomain, contour_style="draw color=white", style="white", xbins=151, ybins=151))
	push!(plots, Plots.Scatter([x[1] for x in xs], [x[2] for x in xs], style="mark=*, mark size=0.5, mark options={draw=white, fill=white}"))
	push!(plots, Plots.Scatter([μ_cma[1]], [μ_cma[2]], style="mark=*, mark size=1.5, mark options={draw=pastelBlue, fill=pastelBlue}"))
	push!(plots, Plots.Scatter([μ_ce[1]], [μ_ce[2]], style="mark=*, mark size=1.5, mark options={draw=pastelRed, fill=pastelRed}"))
	ax = Axis(plots, ymin=ydomain[1], xmin=xdomain[1], ymax=ydomain[2], xmax=xdomain[2], width="5.5cm", height="5.5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end
plot(p)
# -

    using Distributions
    using Random

    p = let
        function flower(x; a=1, b=1, c=4)
			return a*norm(x) + b*sin(c*atan(x[2], x[1]))
		end
        f = x -> flower(x)
        xdomain = ( -3, 3)
		ydomain = ( -3, 3)

        G = GroupPlot(4,4,groupStyle="xlabels at=edge bottom, ylabels at =edge left, horizontal sep=0.5cm, vertical sep=0.5cm", style="xlabel=\$x_1\$, ylabel=\$x_2\$")

        function add_axis!(samples, P)
            plots = Plots.Plot[]
            push!(plots, Plots.Image((x,y)->f([x,y]), xdomain, ydomain, colormap=pasteljet, colorbar=false, xbins=600, ybins=600))
            push!(plots, Plots.Contour(x->pdf(P, x), xdomain, ydomain, contour_style="draw color=white", style="white", xbins=151, ybins=151))
            push!(plots, Plots.Scatter([x[1] for x in samples], [x[2] for x in samples], style="mark=*, mark size=0.5, mark options={draw=white, fill=white}"))
            ax = Axis(plots, ymin=ydomain[1], xmin=xdomain[1], ymax=ydomain[2], xmax=xdomain[2], width="5cm", height="5cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
            push!(G, ax)
        end

        Random.seed!(0)
        σ = 1.0
        x = [2.0, 2.0] # initial point
        μ = copy(x)

        m = 4 + floor(Int, 3*log(length(μ)))
        m_elite = div(m,2)
        cm = 1.0

        n = length(x)
        w′s = log((m+1)/2) .- log.(1:m)
        w′s ./= sum(w′s[1:m_elite])

        μ_eff = 1 / sum(w′s[1:m_elite].^2)

        @assert 1 ≤ μ_eff ≤ m_elite

        cσ = (μ_eff + 2)/(n + μ_eff + 5)
        dσ = 1 + 2max(0, sqrt((μ_eff-1)/(n+1))-1) + cσ
        cΣ = (4 + μ_eff/n)/(n + 4 + 2μ_eff/n)
        c1 = 2/((n+1.3)^2 + μ_eff)
        cμ = min(1-c1, 2*(μ_eff-2+1/μ_eff)/((n+2)^2 + μ_eff))

        α1 = 1 + c1 / cμ
        α2 = 1 + 2μ_eff / (μ_eff + 2)
        α3 = (1 - c1 - cμ)/(n*cμ)
        α_low = min(α1, α2, α3)

        w′_pos =  sum(max(w′, 0) for w′ in w′s)
        w′_neg = -sum(min(w′, 0) for w′ in w′s)
        ws = [w′ ≥ 0 ? w′/w′_pos : α_low*w′/w′_neg for w′ in w′s]

        E_norm = n^0.5*(1-1/(4n)+1/(21*n^2))

        pσ = zeros(n)
        pΣ = zeros(n)
        Σ = Matrix(1.0I, n, n)
        μ = copy(x)

        @assert isapprox(sum(ws[1:m_elite]), 1.0, atol=1e-6)

        for k in 1 : prod(G.dimensions)
            P = MvNormal(μ, σ^2*Σ)
            xs = [rand(P) for i in 1 : m]
            ys = [f(x) for x in xs]
            is = sortperm(ys) # best to worst

            add_axis!(xs, P)

            # selection and mean update
            ys = [(x - μ)/σ for x in xs]
            yw = sum(ws[i]*ys[is[i]] for i in 1 : m_elite)
            μ += cm*σ*yw

            # step-size control
            C = Σ^-0.5
            pσ = (1-cσ)*pσ + sqrt(cσ*(2-cσ)*μ_eff)*C*yw
            σ *= exp(cσ/dσ * (norm(pσ)/E_norm - 1))

            # covariance adaptation
            hσ   = norm(pσ)/sqrt(1-(1-cσ)^(2k)) < (1.4 + 2/(n+1))*E_norm ? 1 : 0
            pΣ = (1-cΣ)*pΣ + hσ*sqrt(cΣ*(2-cΣ)*μ_eff)*yw
            w0 = [ws[i] ≥ 0 ? ws[i] : n*ws[i]/norm(C*ys[is[i]])^2 for i in 1 : m]
            Σ = (1-c1-cμ) * Σ + # regard old matrix # Eq. 47
                c1 * (pΣ*pΣ' + # plus rank one update
                      (1-hσ) * cΣ*(2-cΣ) * Σ) + # minor correction
                cμ * # plus rank mu update
                    sum(w0[i]*ys[is[i]]*ys[is[i]]' for i in 1 : m)

            Σ = triu(Σ)+triu(Σ,1)' # enforce symmetry
        end

        G
    end

    plot(p)

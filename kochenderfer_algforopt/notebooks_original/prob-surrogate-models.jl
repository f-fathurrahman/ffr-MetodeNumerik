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

# # Prob Surrogate Models
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
include("gp.jl")

using LinearAlgebra

p = let

	domx = (-4,4)
	domy = (-4,4)

	g = GroupPlot(3,1,groupStyle="horizontal sep=0.5cm, xlabels at=edge bottom, xticklabels at=edge bottom, ylabels at=edge left, yticklabels at=edge left", style="width=6.5cm, height=6.5cm, axis equal, view={0}{90}, xlabel=\$x_1\$, ylabel=\$x_2\$,")

	Ns = [MvNormal(Float64[0,0], Matrix(1.0I, 2, 2)),
	MvNormal(Float64[0,0], Float64[3 0; 0 0.5]),
	MvNormal(Float64[0,0], [1 0.9; 0.9 1])]
	titles = [L"\mat \Sigma = \begin{bmatrix}1 & 0 \\ 0 & 1\end{bmatrix}",
	L"\mat \Sigma = \begin{bmatrix}3 & 0 \\ 0 & 0.5\end{bmatrix}",
	L"\mat \Sigma = \begin{bmatrix}1 & 0.9 \\ 0.9 & 1\end{bmatrix}"
	]

	for (N,title) in zip(Ns, titles)
	ax = Axis(Plots.Contour((x,y)->pdf(N, [x,y]), domx, domy, levels=[0.001,0.01,0.05,0.1,0.2,0.3], xbins=201, ybins=201),
	xmin=-3, xmax=3, ymin=0, ymax=0.6, title = title, style="contour/labels=false",
	)
	push!(g, ax)
	end

	g
end

plot(p)

# +
using Random

p = let
	g = GroupPlot(3,1,groupStyle="horizontal sep=0.5cm", style="cycle list name = pastelcolors, width=6.5cm")

	Random.seed!(0)
	X = [[x] for x in range(0, stop=10, length=201)]
	function get_plot(m, l, title)
		GP = GaussianProcess(m=m, k=(x, x′) -> exp(-(norm(x - x′)^2/(2*l^2))))

		p = Plots.Plot[]

		for i in 1:5
		    y = rand(GP, X)
		    push!(p, Plots.Linear([x[1] for x in X], y))
		end
		ax = Axis(p, xlabel=L"x",
		             xmin=0, xmax=10, ymin=-3.5, ymax=3.5, title=title)
		if length(g.axes) == 0
		    ax.ylabel = L"y"
		else
		    ax.style="yticklabels={,}"
		end
		ax
	end

	push!(g, get_plot(x -> 0.0, 0.5, L"\ell = \sfrac{1}{2}"))
	push!(g, get_plot(x -> 0.0, 1, L"\ell = 1"))
	push!(g, get_plot(x -> 0.0, 2, L"\ell = 2"))

	g
end

plot(p)

# +
using Random
import SpecialFunctions: gamma, besselk

matern(ν,ℓ) = (x,x′) -> begin
    r = norm(x-x′)
    if isapprox(r, 0.0)
        return 1.0
    end
    g = sqrt(2ν)*r/ℓ
    (2^(1-ν))/gamma(ν) * g^ν * besselk(ν,g)
end

get_nn_kernel(Σ) = (x,x′) -> begin
    x = vcat(1.0, x)
    x′ = vcat(1.0, x′)
    asin((2*(x'*Σ*x′)[1])/sqrt((1+(2*x'*Σ*x)[1])*(1+(2*x′'*Σ*x′)[1])))
end

p = begin
	kernels = Function[]
	titles = LaTeXString[]

	push!(titles, L"Constant: $\sigma_0^2$"); push!(kernels, (x,x′) -> 1.0)
	push!(titles, L"Linear: $\displaystyle \sum_{d=1}^{n} \sigma^2_d x_d  x^\prime_d$"); push!(kernels, (x,x′) -> x⋅x′)
	push!(titles, L"Polynomial: $(x^\top x^\prime + \sigma_0^2)^p$"); push!(kernels, (x,x′) -> (x⋅x′+1)^2)
	push!(titles, L"Exponential: $\exp(-\frac{r}{\ell})$"); push!(kernels, (x,x′) -> exp(-abs(x[1]-x′[1])))
	push!(titles, L"$\gamma$-Exponential: $\exp(-(\frac{r}{\ell})^\gamma)$"); push!(kernels, (x,x′) -> exp(-(norm(x-x′)^0.5)))
	push!(titles, L"Squared Exponential: $\exp(-\frac{r^2}{2\ell^2})$"); push!(kernels, (x,x′) -> exp(-norm(x-x′)^2/2))
	push!(titles, L"Mat\'{e}rn: $\frac{2^{1-\nu}}{\Gamma(\nu)} \left(\sqrt{2\nu}\frac{r}{\ell}\right)^\nu K_\nu \left(\sqrt{2\nu}\frac{r}{\ell}\right)$"); push!(kernels, matern(0.5,1.0))
	push!(titles, L"Rational Quadratic: $\left(1 + \frac{r^2}{2\alpha\ell^2}\right)^{-\alpha}$"); push!(kernels, (x,x′) -> (1+norm(x-x′)^2)^-0.5)
	push!(titles, L"Neural Network: $\sin^{-1}\left(\frac{2\bar{\vect{x}}^\top \Sigma\bar{\vect{x}}^\prime}{\sqrt{(1+2\bar{\vect{x}}^\top\mat \Sigma\bar{\vect{x}})(1+2\bar{\vect{x}}^{\prime\top}\mat\Sigma\bar{\vect{x}}^\prime)}}\right)$"); push!(kernels, get_nn_kernel(Matrix(1.0I, 2, 2)))

	g = GroupPlot(3,3,groupStyle="horizontal sep=0.75cm, vertical sep=1.25cm, xlabels at=edge bottom, xticklabels at=edge bottom, ylabels at=edge left, yticklabels at=edge left",
	                    style="cycle list name = pastelcolors, height=5cm, width=6cm, xmin=-5, xmax=5, ymin=-3.5, ymax=3.5, xlabel=\$x\$, ylabel=\$y\$, every axis title/.style={font=\\footnotesize, at={(0.5,1.15)}, align=center}")

	Random.seed!(0)
	X = [[x] for x in range(-5, stop=5, length=101)]

	function get_plot(title, k)
	    p = Plots.Plot[]
	    GP = GaussianProcess(k=k)
	    for i = 1:5
	        y = rand(GP, X)
	        push!(p, Plots.Linear([x[1] for x in X], y))
	    end
	    Axis(p, title=title)
	end

	for (title, k) in zip(titles, kernels)
	    push!(g, get_plot(title, k))
	end

	g
end

plot(p)

# +
using Random

p = begin
	g = GroupPlot(3,1,groupStyle="horizontal sep=0.75cm, vertical sep=1.5cm, xlabels at=edge bottom, ylabels at=edge left",
	                       style="width=6cm, height=6cm, xlabel=\$x_1\$, ylabel=\$x_2\$, xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")

	x_arr = range(0, stop=10, length=71)
	y_arr = range(0, stop=10, length=71)

	for l in [0.5, 1.0, 2.0]

		GP = GaussianProcess(k = (x, x′) -> exp(-(norm(x - x′)^2/(2*l^2))))

	    # create random "true" function
	    Random.seed!(0)
	    X = vec([[x, y] for x in x_arr, y in y_arr])
	    z = rand(GP, X)
	    n = length(x_arr)
	    A = Float64[z[(i-1)*n+j] for i in 1:length(x_arr), j in 1 : length(y_arr)]

	    push!(g, Axis(Plots.Contour(A, x_arr, y_arr)))
	end

	g.axes[1].title = L"\ell = \sfrac{1}{2}"
	g.axes[2].title = L"\ell = 1"
	g.axes[3].title = L"\ell = 2"

	g
end

plot(p)

# +
using Random

p = let

	# create random "true" function
	Random.seed!(0)
	GP = GaussianProcess()
	x_arr = range(0, stop=10, length=101)
	X = [[x] for x in collect(x_arr)]
	y = rand(GP, X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true objective function")

	# create training set
	for i in [5, 20, 35, 45]
		push!(GP, X[i], y[i])
	end
	p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="solid, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	# predict
	μₚ, νₚ = predict(GP, X)
	p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted function mean")
	σₚ = sqrt.(νₚ)
	upperConfidence = μₚ + 1.96*σₚ
	lowerConfidence = μₚ - 1.96*σₚ
	p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
	p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
	p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence interval}")

	Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_pred_μ, p_train],
		xmin=0, xmax=8, style="axis on top, ytick=\\empty, xtick=\\empty",
		xlabel=L"x", ylabel=L"y", width="7.5cm", legendPos="outer north east",
	)
end

plot(p)
# -

        struct GaussianProcessWithGradient
            kff::Function
            k∇f::Function
            k∇∇::Function
            X::Vector{Vector{Float64}}
            y::Vector{Float64}
            ∇y::Vector{Vector{Float64}}
        end
        function GaussianProcessWithGradient(α::Float64, γ::Float64)

            squaregammanorm = (x, γ) -> (x.*γ)⋅x

            kff = (x,x′) -> α*exp(-0.5*squaregammanorm(x-x′, γ))
            k∇f = (x,x′,j) -> -γ*α*(x′[j]-x[j])*exp(-0.5*squaregammanorm(x-x′, γ))
            k∇∇ = (x,x′,i,j) -> γ*α*((i == j ? 1 : 0) - γ*(x[j]-x′[j])*(x[i]-x′[i]))*exp(-0.5*squaregammanorm(x-x′, γ))

            return GaussianProcessWithGradient(
                kff,
                k∇f,
                k∇∇,
                Vector{Float64}[],
                Float64[],
                Vector{Float64}[],
            )
        end

        function Base.push!(GP::GaussianProcessWithGradient, x::Vector{Float64}, y::Real, ∇y::Vector{Float64})
            push!(GP.X, x)
            push!(GP.y, y)
            push!(GP.∇y, ∇y)
            return GP
        end
        function Base.pop!(GP::GaussianProcessWithGradient)
            pop!(GP.X)
            pop!(GP.y)
            pop!(GP.∇y)
            return GP
        end

        K(X::Dset, X′::Dset, k::Function) = [k(x,x′) for x in X, x′ in X′]
        K(X::Dset, X′::Dset, k::Function, i::Int) = [k(x,x′,i) for x in X, x′ in X′]
        K(X::Dset, X′::Dset, k::Function, i::Int, j::Int) = [k(x,x′,i,j) for x in X, x′ in X′]

        function K∇f(X::Dset, X′::Dset, k::Function)
            n = length(X)
            m = length(X′)
            d = length(X[1])
            retval = Array{Float64}(undef, n, m*d)
            for (i,x) in enumerate(X)
                for (j,x′) in enumerate(X′)
                    for a in 1 : d
                        index = d*(j-1) + a
                        retval[i, index] = k(x,x′,a)
                    end
                end
            end
            return retval
        end
        function K∇∇(X::Dset, X′::Dset, k::Function)
            n = length(X)
            m = length(X′)
            d = length(X[1])
            retval = Array{Float64}(undef, n*d, m*d)
            for (i,x) in enumerate(X)
                for (j,x′) in enumerate(X′)
                    for a in 1 : d
                        for b in 1 : d
                            retval[d*(i-1) + a, d*(j-1) + b] = k(x,x′,a,b)
                        end
                    end
                end
            end
            return retval
        end

        function predict(GPN::GaussianProcessWithGradient, X_pred::Dset)
            _Kff = K(GPN.X, GPN.X, GPN.kff)
            _K∇f = K∇f(GPN.X, GPN.X, GPN.k∇f)
            _K∇∇ = K∇∇(GPN.X, GPN.X, GPN.k∇∇)

            _Kpf = K(X_pred, GPN.X, GPN.kff)
            _K∇p = K∇f(GPN.X, X_pred, GPN.k∇f)
            _Kpp = K(X_pred, X_pred, GPN.kff)

            v = hcat(_Kpf, _K∇p')'
            V = vcat(hcat(_Kff, _K∇f'), hcat(_K∇f, _K∇∇))

            μₚ = v'/V * vcat(GPN.y, -vcat(GPN.∇y...))
            νₚ = diag(_Kpp - v'/V * v) .+ eps()
            return (μₚ, νₚ)
        end

        p = let
            α = 1.0
            γ = 1.5

            domain = (0,8)
            f = x -> 1 - exp(-(x[1]-4)^2)
			∇f = x -> [2*exp(-(x[1]-4)^2)*(x[1]-4)]

            squaregammanorm = (x, γ) -> (x.*γ)⋅x

            GP = GaussianProcess(k=(x,x′) -> α*exp(-0.5*squaregammanorm(x-x′,γ)))
            GPN = GaussianProcessWithGradient(α, γ)
            x_arr = collect(range(domain[1], stop=domain[2], length=151))
            X = [[x] for x in collect(x_arr)]
            for x in [[1.0], [2.0], [3.5], [4.5]]
                push!(GP, x, f(x))
                push!(GPN, x, f(x), ∇f(x))
            end

            p_train = Plots.Scatter([x[1] for x in GPN.X], GPN.y, style="draw opacity=0, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")
            p_true = Plots.Linear(x_arr, f.(X), style="solid, black, mark=none", legendentry="true function")

            # predict w/o gradient
            μₚ, νₚ = predict(GP, X)
            p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted mean without gradient")
            σₚ = sqrt.(νₚ)
            upperConfidence = μₚ + 1.96*σₚ
            lowerConfidence = μₚ - 1.96*σₚ
            p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
            p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
            p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence interval without gradient}")

            # predict with gradient
            ∇μₚ, ∇νₚ = predict(GPN, X)
            p_pred_∇μ = Plots.Linear(x_arr, ∇μₚ, style="solid, thick, pastelRed, mark=none", legendentry="predicted mean with gradient")
            ∇σₚ = sqrt.(∇νₚ)
            upperConfidence∇ = ∇μₚ + 1.96*∇σₚ
            lowerConfidence∇ = ∇μₚ - 1.96*∇σₚ
            p_pred_∇hi = Plots.Linear(x_arr, upperConfidence∇, style="draw=none, mark=none, name path=A, forget plot")
            p_pred_∇lo = Plots.Linear(x_arr, lowerConfidence∇, style="draw=none, mark=none, name path=B, forget plot")
            p_pred∇ = Plots.Command("\\addplot[pastelRed!70] fill between[of=A and B];\\addlegendentry{confidence interval with gradient}")

            Axis([p_pred_hi, p_pred_lo, p_pred, p_pred_∇hi, p_pred_∇lo, p_pred∇, p_pred_μ, p_pred_∇μ, p_train, p_true],
              style="enlarge x limits=0, axis on top, ytick=\\empty, xtick=\\empty",
              xlabel=L"x", ylabel=L"y", xmin=0, xmax=8, width="7.5cm", legendPos="outer north east"
            )
        end
        plot(p)

# +
using Random

p = let

	GP = GaussianProcess(
		m = x -> 0.0,
		k = (x, x′) -> exp(-(norm(x - x′))^2), # squared exponential
		ν = 0.05,
	)

	# create random "true" function
	Random.seed!(0)
	x_arr = range(0, stop=10, length=101)
	X = [[x] for x in collect(x_arr)]
	y = rand(GP, X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true objective function")

	# create training set
	for i in [5, 20, 35, 45, 45]
		push!(GP, X[i], y[i] + randn().*sqrt(GP.ν))
	end
	p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="draw opacity=0, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	# predict
	μₚ, νₚ = predict(GP, X)
	σₚ = sqrt.(νₚ)
	p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted function mean")
	upperConfidence = μₚ + 1.96*σₚ
	lowerConfidence = μₚ - 1.96*σₚ
	p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
	p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
	p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

	Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_pred_μ, p_train],
		xmin=0, xmax=8, style="axis on top, ytick=\\empty, xtick=\\empty",
		xlabel=L"x", ylabel=L"y", width="7.5cm", legendPos="outer north east",
	)
end

plot(p)
# -



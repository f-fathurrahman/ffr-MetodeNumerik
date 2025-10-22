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

# # Uncertaintyprop
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

using LinearAlgebra
using Distributions
using Optim
using ForwardDiff
import QuadGK: quadgk
p = let
	f = (x,z1,z2) -> sin(x+z1)*cos(x+z2)
	dom =  (-π,π)

	function taylor_approx(f, μ, ν, secondorder=false)
		μhat = f(μ)
		∇ = (z -> ForwardDiff.gradient(f, z))(μ)
		νhat = ∇.^2⋅ν
		if secondorder
			H = (z -> ForwardDiff.hessian(f, z))(μ)
			μhat += (diag(H)⋅ν)/2
			νhat += ν⋅(H.^2*ν)/2
		end
		return (μhat, νhat)
	end

	μ_ = [0.0,0.0]
    ν = [0.1,0.2]
    P1 = Normal(μ_[1],sqrt(ν[1]))
    P2 = Normal(μ_[2],sqrt(ν[2]))
    P = MvNormal(μ_, sqrt.(ν))

    ϵ = 1e-6
    E = x -> quadgk(z1->quadgk(z2->f(x,z1,z2)  *pdf(P2,z2), -Inf, Inf, atol=ϵ)[1]*pdf(P1,z1), -Inf,Inf, atol=ϵ)[1]
    V = x -> quadgk(z1->quadgk(z2->f(x,z1,z2)^2*pdf(P2,z2), -Inf, Inf, atol=ϵ)[1]*pdf(P1,z1), -Inf,Inf, atol=ϵ)[1] - E(x)^2

    G = GroupPlot(1,2, groupStyle="vertical sep=0.25cm, xlabels at=edge bottom, xticklabels at=edge bottom", style="width=10cm, height=3cm, xmin=$(convert(Float64, dom[1])), xmax=$(convert(Float64, dom[2])), xlabel=\$x\$")
    push!(G,
    Axis([
        Plots.Linear(E, dom, style="solid, thick, black, mark=none"),
        Plots.Linear(x->taylor_approx(z->f(x,z[1],z[2]), μ_, ν)[1], dom, style="solid, pastelBlue, mark=none"),
        Plots.Linear(x->taylor_approx(z->f(x,z[1],z[2]), μ_, ν, true)[1], dom, style="solid, pastelRed, mark=none"),
        ], ylabel="mean"))
    push!(G,
    Axis([
        Plots.Linear(V, dom, style="solid, thick, black, mark=none", legendentry="true"),
        Plots.Linear(x->taylor_approx(z->f(x,z[1],z[2]), μ_, ν)[2], dom, style="solid, pastelBlue, mark=none", legendentry="first-order Taylor approx."),
        Plots.Linear(x->taylor_approx(z->f(x,z[1],z[2]), μ_, ν, true)[2], dom, style="solid, pastelRed, mark=none", legendentry="second-order Taylor approx."),
        ], ylabel="variance", ymin=0, style="scaled ticks=false, tick label style={/pgf/number format/fixed},legend style={draw=none, at={(0.5,-0.6)}, anchor=north, legend columns=-1, /tikz/every even column/.append style={column sep=0.5cm}}"))
end
plot(p)

using Optim
using Distributions
p = let

	f = (x,z1,z2) -> sin(x+z1)*cos(x+z2) + 2
	dom =  (-0.9π,0.9π)

	ν1 = 1.0
	ν2 = 0.5
	P1 = Normal(0,sqrt(ν1))
	P2 = Normal(0,sqrt(ν2))

	n = 2
	k = 1
	ϕ = (x,z1,z2) ->  -k*log(f(x,z1,z2)) - logpdf(P1, z1) - logpdf(P2,z2)

	get_z_star = x -> optimize(z->ϕ(x,z[1],z[2]), zeros(n), Optim.BFGS()).minimizer

	∂1 = (x,z1,z2) ->  -k*cos(x+z1)*cos(x+z2) / f(x,z1,z2) + z1
	∂2 = (x,z1,z2) ->   k*sin(x+z1)*sin(x+z2) / f(x,z1,z2) + 2z2
	∂11 = (x,z1,z2) -> k*cos(x+z1)^2*cos(x+z2)^2 / f(x,z1,z2)^2 + k*sin(x+z1)*cos(x+z2)/f(x,z1,z2) + 1
	∂22 = (x,z1,z2) -> k*sin(x+z1)^2*sin(x+z2)^2 / f(x,z1,z2)^2 + k*sin(x+z1)*cos(x+z2)/f(x,z1,z2) + 2
	∂12 = (x,z1,z2) -> k*cos(x+z1)*sin(x+z2)/f(x,z1,z2) - k*sin(x+z1)*cos(x+z1)*sin(x+z2)*cos(x+z2)/f(x,z1,z2)^2

	H = (x,z1,z2) -> [∂11(x,z1,z2) ∂12(x,z1,z2); ∂12(x,z1,z2) ∂22(x,z1,z2)]

	ϵ = 1e-5
	E = x -> quadgk(z1->quadgk(z2->f(x,z1,z2)  *pdf(P2,z2), -Inf, Inf, atol=ϵ)[1]*pdf(P1,z1), -Inf,Inf, atol=ϵ)[1]
	_E = x -> begin
	    zstar = get_z_star(x)
	    Hstar = H(x, zstar[1], zstar[2])
	    dHstar = det(Hstar)
	    if dHstar > 0.0
	        return clamp(exp(-ϕ(x, zstar[1], zstar[2]))*(2π)^(n/2)*dHstar^(-0.5), 1.0, 10.0)
	    else
	        return 0.0
	    end
	end

	Axis(
	    [
	    Plots.Linear(E, dom, style="solid, black, mark=none", xbins=101, legendentry="true")
	    Plots.Linear(_E, dom, style="solid, pastelBlue, mark=none", xbins=251, legendentry="Laplace approx.")
	     ], xlabel=L"x", ylabel=L"\Expectation[x]", width="8cm", height="6cm", style="enlarge x limits=false, legend cell align=left, legend style={draw=none}, legend pos=outer north east",
	)
end
plot(p)

# +
using Distributions
using Polynomials
using Random
import QuadGK: quadgk

p = let
	function hermite(i)
	    p = Polynomial([1])
	    x = Polynomial([0,1])
	    for j in 2 : i
	        p = x*p - derivative(p)
	    end
	    return p
	end

	bs = [hermite(i) for i in 1 : 4]

	Z = Normal(0.0, 1.0)

	function fobj(x::Float64, z::Float64=0.0)
	    return 1 - exp(-(x+z-1)^2) - 2exp(-(x+z-3)^2)
	end
	E = x -> quadgk(z->fobj(x,z)*pdf(Z,z), -Inf,Inf)[1]

	function polyE(x::Float64; m::Int=7)
	    zs = rand(Z,m)
	    ys = [fobj(x,z) for z in zs]
	    Bmat = [b(z) for z in zs, b in bs]
	    θ = pinv(Bmat)*ys
	    # return the mean
	    return θ[1]
	end

	function polyEBounds(x::Float64, m::Int=7, m2::Int=1000)
	    N = fit_mle(Normal, [polyE(x, m=m) for i in 1 : m2])
	    return (quantile(N, 0.05), mean(N), quantile(N, 0.95))
	end

	dom = (-2,6)

	sample_counts = [10,30,50]
	G = GroupPlot(1, length(sample_counts), groupStyle="vertical sep=0.5cm, xlabels at=edge bottom, xticklabels at=edge bottom",
						style="width=9.5cm, height=3cm, xlabel=\$x\$, ylabel={\$\\Expectation[f \\mid x]\$}, enlarge x limits=0, ymax=1.5, ymin=-1.5")

	x_arr = range(dom[1], stop=dom[2], length=200)
	for (i,m) in enumerate(sample_counts)
		Random.seed!(0)
		p = Plots.Plot[]
		push!(p, Plots.Linear(fobj, dom, xbins=100, style="solid, pastelRed, mark=none", legendentry= i==1 ? "noise-free" : nothing))
		push!(p, Plots.Linear(E, dom, xbins=100, style="solid, black, mark=none", legendentry= i==1 ? "exact" : nothing))

		bounds = [polyEBounds(x, m, 1000) for x in x_arr]
		push!(p, Plots.Linear(x_arr, [b[1] for b in bounds], style="name path=A, draw=none, mark=none, forget plot"))
		push!(p, Plots.Linear(x_arr, [b[3] for b in bounds], style="name path=B, draw=none, mark=none, forget plot"))
		push!(p, Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];"))
		if i == 1
			push!(p, Plots.Command("\\addlegendentry{\\SI{95}{\\percent} interval}"))
		end
		push!(p, Plots.Linear(x_arr, [b[2] for b in bounds], style="solid, pastelBlue, mark=none", legendentry= i==1 ? "mean" : nothing))
		push!(p, Plots.Node("\\scriptsize \$$m\$ samples", 0.02, 0.25, style="right", axis="axis description cs"))
		push!(G, Axis(p))
	end

	G.axes[1].style = "legend pos=outer north east"

	G
end
plot(p)

# +
using Polynomials
function legendre(i)
	n = i-1
    p = Polynomial([-1,0,1])^n
    for i in 1 : n
        p = derivative(p)
    end
    return p / (2^n * factorial(n))
end
function laguerre(i)
    p = Polynomial([1])
    for j in 2 : i
        p = integrate(derivative(p) - p)
    end
    return p
end
function hermite(i)
    p = Polynomial([1])
    x = Polynomial([0,1])
    for j in 2 : i
        p = x*p - derivative(p)
    end
    return p
end

p = begin
	G = GroupPlot(3,1,groupStyle="horizontal sep=0.75cm, vertical sep=1.5cm, xlabels at=edge bottom, xticklabels at=edge bottom",
	                       style="cycle list name = pastelcolors, width=6cm, every axis title/.style={font=\\footnotesize, at={(0.5,1.1)}, align=center}")

	p = Plots.Plot[]
	for i in 1 : 6
		L = legendre(i)
	    push!(p, Plots.Linear(x->L(x), (-1, 1), xbins=201, style="solid, mark=none"))
	end
	push!(G, Axis(p, xlabel=L"z", title="Legendre", xmin=-1, xmax=1))

	p = Plots.Plot[]
	for i in 1 : 6
		L = laguerre(i)
	    push!(p, Plots.Linear(x->L(x), (0, 15), xbins=201, style="solid, mark=none"))
	end
	push!(G, Axis(p, xlabel=L"z", title="Laguerre", ymin=-20, ymax=20, xmin=0, xmax=15))

	p = Plots.Plot[]
	for i in 1 : 6
	    H = hermite(i)
	    push!(p, Plots.Linear(z->clamp(H(z), -17, 17), (-5, 5), xbins=201, style="solid, mark=none", legendentry="\$b_$i\$"))
	end
	push!(G, Axis(p, xlabel=L"z", style="legend pos=outer north east", title="Hermite", xmin=-5, xmax=5, ymin=-15, ymax=15))

	G
end

plot(p)

# +
using QuadGK
import QuadGK: quadgk
import Printf: @sprintf

p = let
	f = z -> sin(z*π)
	X = [[-1], [-0.2], [0.3], [0.7], [0.9]]
	y = [f(x[1]) for x in X]

	var_term(b, ρ, dom) = quadgk(x->b(x)^2*ρ(x), dom...)[1]
	var_term_legendre(i) = var_term(legendre(i), x->1/2, [-1,1])
	function get_var(θs, var_terms)
	    μ_ = θs[1]
	    return θs.^2⋅var_terms - μ_^2
	end


	function get_plot(i)
	    bases = [z->legendre(j)(z[1]) for j in 1 : i]
	    var_terms = [var_term_legendre(j) for j in 1 : i]
	    B = [b(x) for x in X, b in bases]
	    θ = pinv(B)*y
	    fhat = z -> sum(θ[j] * bases[j](z) for j in 1 : i)
	    ν = get_var(θ, var_terms)
	    return Plots.Linear(z->fhat([z]), (-1,1), style="solid, mark=none", legendentry="\$i = $i \\quad \\mu = $(@sprintf("%+6.3f", θ[1])) \\quad \\nu = $(@sprintf("%6.3f", ν))\$")
	end

	p = Plots.Plot[]
	push!(p, Plots.Linear(f, (-1,1), style="solid, black, mark=none", legendentry="true"))
	for i in 1 : 5
	    push!(p, get_plot(i))
	end
	push!(p, Plots.Scatter([x[1] for x in X], y, style="solid, mark=*, mark size=1, mark options={draw=black, fill=black}"))
	Axis(p, style="legend pos=outer north east, cycle list name = pastelcolors", xlabel=L"z", width="7.5cm", xmin=-1, xmax=1)
end
plot(p)

# +
using Distributions
using Polynomials

p = let
	f = z -> sin(z*π)
	X = [[2.1], [2.5], [3.3], [4.7], [3.9]]
	y = [f(x[1]) for x in X]

	dom = (2,5)
	P = TruncatedNormal(3, 1, dom...)
	ρ = z -> pdf(P, z)

	var_term(b, ρ, dom) = begin
		b2 = b*b
		quadgk(x->b2(x)*ρ(x), dom...)[1]
	end
	function get_var(θs, var_terms)
	    μ_ = θs[1]
	    return θs.^2⋅var_terms - μ_^2
	end
	function orthogonal_recurrence(bs, p, dom, ϵ=1e-6)
	    i = length(bs)
	    c1 = quadgk(z->z*bs[i](z)^2*p(z), dom..., atol=ϵ)[1]
	    c2 = quadgk(z->  bs[i](z)^2*p(z), dom..., atol=ϵ)[1]
	    α = c1 / c2
	    if i > 1
	    	c3 = quadgk(z->bs[i-1](z)^2*p(z), dom..., atol=ϵ)[1]
	    	β = c2 / c3
	    	return Polynomial([-α, 1])*bs[i] - β*bs[i-1]
	    else
	    	return Polynomial([-α, 1])*bs[i]
		end
	end

	bs = [Polynomial([1.0])]
	var_terms = [var_term(bs[end], ρ, dom)]

	p = Plots.Plot[]
	push!(p, Plots.Linear(f, dom, style="solid, black, mark=none", legendentry="true"))
	for i in 1 : 5

	    B = [b(x[1]) for x in X, b in bs]
	    θ = pinv(B)*y
	    fhat = z -> sum(θ[j] * bs[j](z) for j in 1 : i)
	    ν = get_var(θ, var_terms)

	    push!(p, Plots.Linear(z->fhat(z), dom, style="solid, mark=none", legendentry="\$i = $i \\quad \\mu = $(@sprintf("%+6.3f", θ[1])) \\quad \\nu = $(@sprintf("%6.3f", ν))\$"))

	    push!(bs, orthogonal_recurrence(bs, ρ, dom))
	    push!(var_terms, var_term(bs[end], ρ, dom))
	end
	push!(p, Plots.Scatter([x[1] for x in X], y, style="solid, mark=*, mark size=1, mark options={draw=black, fill=black}"))
	Axis(p, style="legend pos=outer north east, cycle list name = pastelcolors", xlabel=L"z", width="7.5cm", xmin=dom[1], xmax=dom[2])
end
plot(p)

# +
if !@isdefined(GaussianProcess)
	include("../chapter/gp.jl")
end

p = let
	f = (x,z) -> sin(x+z[1])*cos(x+z[2])
	w = [1.0,1.0]
	μz = [0.0,0.0]
	Σz = Matrix(Diagonal([1.0,0.5]))
	D = MvNormal(μz, Σz)

	function bayesian_monte_carlo(GP, w, μz, Σz)
		W = Matrix(Diagonal((w.^2)))
		invK = inv(K(GP.X, GP.X, GP.k))
		q = [exp(-0.5*((z-μz)'inv(W+Σz)*(z-μz))[1]) for z in GP.X]
		q .*= (det(W\Σz + I))^(-0.5)
		μ = q'*invK*GP.y
		ν = (det(2W\Σz + I))^(-0.5) - (q'*invK*q)[1]
		return (μ, ν)
	end

	dom =  (-0.9π,0.9π)
	ϵ = 1e-5
	P1 = Normal(0,sqrt(1.0))
	P2 = Normal(0,sqrt(1/2))
	E = x -> QuadGK.quadgk(z1->QuadGK.quadgk(z2->f(x,[z1,z2])  *pdf(P2,z2), -Inf, Inf, atol=ϵ)[1]*pdf(P1,z1), -Inf,Inf, atol=ϵ)[1]

	x_arr = collect(range(dom[1], stop=dom[2], length=251))
	samples = [[rand(D) for j in 1 : 10] for i in 1 : length(x_arr)]

	y_bmc = Array{Float64}(undef, length(x_arr))
	for i in 1 : length(x_arr)
		x = x_arr[i]
		GP = GaussianProcess(k = (x,x′) ->
		    exp(-0.5*sum((x[i]-x′[i])^2/w[i]^2 for i in 1 : length(x)))
		)
		for z in samples[i]
		    push!(GP, z, f(x,z))
		end
		y_bmc[i] = bayesian_monte_carlo(GP, w, μz, Σz)[1]
	end

	y_mean = [mean(f(x,z) for z in zs) for (x,zs) in zip(x_arr, samples)]

	Axis(
	    [
	    Plots.Linear(E, dom, style="solid, black, mark=none", xbins=101, legendentry="true")
	    Plots.Linear(x_arr, y_mean, style="solid, pastelRed, mark=none", legendentry="sample mean")
	    Plots.Linear(x_arr, y_bmc, style="solid, pastelBlue, mark=none", legendentry="Bayesian MC")
	     ], xlabel=L"x", ylabel=L"\Expectation[f \mid x]", width="8cm", height="6cm", style="enlarge x limits=false, legend cell align=left, legend style={draw=none}, legend pos=outer north east",
	)
end
plot(p)

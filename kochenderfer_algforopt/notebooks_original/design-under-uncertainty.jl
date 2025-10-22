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

# # Design Under Uncertainty
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

	using Distributions
	using Optim
	import QuadGK: quadgk

	p = let

		a = 1
		b = 2
		c = 4.5
		f = x -> 2.5 - 2exp(-100(x-a)^2) - exp(-(x-b)^2)  - 1.4exp(-0.15(x-c)^2)

		p = Plots.Plot[]
		x_arr = append!(collect(range(0.0,stop=1.5,length=300)), collect(range(1.5,stop=9,length=200)[2:end]))
		push!(p, Plots.Linear(x_arr, f.(x_arr), style="solid, thick, black, mark=none", legendentry="noise-free"))
		labels=["very low", "low", "high", "very high"]
		for (i,σ) in enumerate([0.1,0.5,1.0,2.0])
		    P = Normal(0,σ)
		    g = x -> quadgk(z->f(x+z)*pdf(P, z), -5, 5)[1]
		    push!(p, Plots.Linear(g, (0, 9), xbins=401, style="solid, pastelBlue!$(100-25(i-1)), mark=none", legendentry=labels[i]*" noise"))
		end

		min_a = optimize(f, a-0.5, a+0.5).minimizer
		min_b = optimize(f, b-0.5, b+0.5).minimizer
		min_c = optimize(f, c-1.0, c+1.0).minimizer

		Axis(p, width="12cm", height="6cm", style="axis lines=left, xtick={$(min_a), $(min_b), $(min_c)}, xticklabels={\$a\$, \$b\$, \$c\$},  xticklabel style={text height=2ex}, ytick=\\empty, legend style={draw=none, at={(0.5,-0.25)}, anchor=north, legend columns=-1, /tikz/every even column/.append style={column sep=0.5cm}}",
		    xlabel=L"x", ylabel="expected value", xmin=-0.5, ymin=-0.25)
    end
    plot(p)

using Optim
p = let
	f = x -> -x*(x < 0) + x^2*(x > 0)
	dom = (-1,1)
	p = Plots.Plot[]
	push!(p, Plots.Linear(f, dom, style="solid, black, mark=none", legendentry="true"))

	F_lin(x,ϵ) = -x + ϵ    # valid for x - ϵ ≤ 0
	F_qua(x,ϵ) = (x + ϵ)^2 # valid for x + ϵ > 0

	for (i,ϵ) in enumerate([0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0])
		x_arr = Float64[dom[1]]
		y_arr = Float64[F_lin(x_arr[1], ϵ)]
		x_vertex = -(2ϵ+1)/2 + 0.5*sqrt((2ϵ+1)^2 - 4*(ϵ^2-ϵ))
		for x in range(x_vertex, stop=dom[2], length=101)
			push!(x_arr, x)
			push!(y_arr, F_qua(x, ϵ))
		end
		push!(p, Plots.Linear(x_arr, y_arr, style="solid, pastelBlue!$(100-9(i-1)), mark=none", legendentry="\$\\epsilon = $ϵ\$"))
	end
	Axis(p, xmin=-1, xmax=1, ymax=4, xlabel=L"x", ylabel=L"y", style="width=8cm, legend pos=outer north east")
end
plot(p)

p = let
	Z = ϵ -> [-ϵ, ϵ]
	f = (x, z) -> (x+z)^2 + 6exp(-(x+z)^2)

	g = (x,ϵ) -> maximum(f(x,z) for z in range(Z(ϵ)[1], stop=Z(ϵ)[2], length=51))
	dom = (-3,3)

	G = GroupPlot(1,2, groupStyle="vertical sep=0.5cm, xlabels at=edge bottom, xticklabels at=edge bottom", style="xlabel=\$x\$, ylabel=\$y\$, xmin=$(dom[1]), xmax=$(dom[2]), ymin=2, ymax=10, width=9.5cm, height=4cm")

	p = Plots.Plot[]
	for (i,ϵ) in enumerate(range(0.0,stop=3.0,length=11))
		if i == 1
			push!(p, Plots.Linear(x->g(x,ϵ), dom, style="solid, black, mark=none", legendentry="noise-free"))
		else
			push!(p, Plots.Linear(x->g(x,ϵ), dom, xbins=101, style="solid, pastelBlue!$(100-9(i-2)), mark=none", legendentry="\$\\epsilon = $(ϵ)\$"))
		end
	end
	push!(p, Plots.Linear([ 2,  2],[2,10], style="solid, gray, mark=none, forget plot"))
	push!(p, Plots.Linear([-2, -2],[2,10], style="solid, gray, mark=none, forget plot"))
	push!(p, Plots.Node(L"\mathcal{Z}",0,8.1))
	push!(p, Plots.Node(L"x^*",0,5))
	push!(p, Plots.Linear([-2, 2],[7.5,7.5], style="solid, thick, red, mark=none, forget plot"))
	push!(p, Plots.Linear([-2,-2],[7.25,7.75], style="solid, thick, red, mark=none, forget plot"))
	push!(p, Plots.Linear([ 2, 2],[7.25,7.75], style="solid, thick, red, mark=none, forget plot"))
	push!(p, Plots.Scatter([0],[f(0,0)], style="solid, black, mark=*, mark size=1.5, mark options={solid, draw=black, fill=black}, forget plot"))
	push!(G, Axis(p, style="legend pos=outer north east, legend style={draw=none, fill=none}"))

	p = Plots.Plot[]
	for (i,ϵ) in enumerate(range(0.0,stop=3.0,length=11))
		if i == 1
			push!(p, Plots.Linear(x->g(x,ϵ), dom, style="solid, black, mark=none"))
		else
			push!(p, Plots.Linear(x->g(x,ϵ), dom, xbins=201, style="solid, pastelBlue!$(100-9(i-2)), mark=none"))
		end
	end
	push!(p, Plots.Linear([ 2,  2],[2,10], style="solid, gray, mark=none"))
	push!(p, Plots.Linear([-2, -2],[2,10], style="solid, gray, mark=none"))
	push!(p, Plots.Linear([dom[1], dom[2]],[5,5], style="solid, gray, mark=none"))

	z_hi = -0.478808
	z_lo = -2
	xstar = (z_hi + z_lo)/2

	push!(p, Plots.Node(L"\mathcal{Z}",-1.25,8))
	push!(p, Plots.Node(L"x^*",xstar+0.5,f(xstar,0)+0.2))
	push!(p, Plots.Linear([z_lo,z_hi],[7.5,7.5], style="solid, thick, red, mark=none"))
	push!(p, Plots.Linear([z_lo,z_lo],[7.25,7.75], style="solid, thick, red, mark=none"))
	push!(p, Plots.Linear([z_hi,z_hi],[7.25,7.75], style="solid, thick, red, mark=none"))
	push!(p, Plots.Scatter([xstar],[f(xstar,0)], style="solid, black, mark=*, mark size=1.5, mark options={solid, draw=black, fill=black}"))

	push!(G, Axis(p))
end
plot(p)

using Distributions
p  = let
	f = x -> sin(2x)/x
	dom = (-10,10)
	p = Plots.Plot[]
	push!(p, Plots.Linear(f, dom, xbins=151, style="solid, black, mark=none", legendentry="\$\\nu = 0\$"))
	for (i,ϵ) in enumerate([0.5,1,1.5,2])
	    P = Normal(0, ϵ)
	    F = x -> quadgk(z->pdf(P, z)*f(x+z), -10, 10)[1]
	    push!(p, Plots.Linear(F, dom, style="smooth,solid, pastelBlue!$(100-10(i-1)), mark=none", legendentry="\$\\nu = $ϵ\$"))
	end
	Axis(p, xmin=dom[1], xmax=dom[2], xlabel=L"x", ylabel=L"y", style="width=9cm, height=4cm, legend cell align=left, legend pos=outer north east")
end
plot(p)

using Distributions
p  = let
	f = x -> isapprox(x, 0.0) ? 2.0 : sin(2x)/x
	dom = (-10,10)
	P = Normal(0, 0.1)
	p = Plots.Plot[]
	push!(p, Plots.Linear(f, dom, xbins=151, style="solid, thick, black, mark=none", legendentry="noise-free"))
	# for (i,k) in enumerate([0.1, 0.5, 1.0, 1.5, 2.0])
	for (i,k) in enumerate([1.0, 1.25, 1.5, 1.75, 2.0])
	    F = x -> quadgk(z->pdf(P, z)*sign(f(x+z))*abs(f(x+z))^k, -Inf, Inf)[1]
	    push!(p, Plots.Linear(F, dom, xbins=150, style="solid, pastelBlue!$(100-15(i-1)), mark=none", legendentry="\$k = $k\$"))
	end
	Axis(p, xmin=dom[1], xmax=dom[2], xlabel=L"x", ylabel=L"y", style="width=11cm, height=6cm, legend cell align=left, legend style={draw=none, at={(0.5,-0.15)}, anchor=north, legend columns=1},")
end
plot(p)

# +
using Distributions
using Optim

p = let

	f = x -> x^3 + 1
	p = Plots.Plot[]
	push!(p, Plots.Linear(f, (0, 5), style="solid, black, mark=none"))

	ymin = -3

	δ = 0.49
	P = Normal(1.5,0.25)
	dom = quantile.(P, [0.5-δ, 0.5+δ])
	a, b = dom[1], dom[2]
	fa, fb = f(a), f(b)
	push!(p, Plots.Linear(x->2pdf(P, x)+ymin, (a, b), style="solid, pastelBlue, mark=none"))
	push!(p, Plots.Linear([-0.5, a, a], [fa, fa, ymin], style="draw=none, mark=none, name path=A, forget plot"))
	push!(p, Plots.Linear([-0.5, b, b], [fb, fb, ymin], style="draw=none, mark=none, name path=B, forget plot"))
	push!(p, Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];"))
	_μ = (fa+fb)/2
	σ = optimize(σ->(quantile(Normal(_μ, σ), 0.5+δ) - fb)^2, 0.1, 10.0).minimizer
	Q = Normal(_μ, σ)
	y_arr = range(fa, stop=fb, length=101)
	x_arr = [3*(pdf(Q, y)-pdf(Q,fa))-0.5 for y in y_arr]
	push!(p, Plots.Linear(x_arr, y_arr, style="solid, pastelBlue, mark=none"))
	push!(p, Plots.Node("robust",0,fb, style="above right"))

	P = Normal(4,0.25)
	dom = quantile.(P, [0.5-δ, 0.5+δ])
	a, b = dom[1], dom[2]
	fa, fb = f(a), f(b)
	push!(p, Plots.Linear(x->2pdf(P, x)+ymin, (a, b), style="solid, pastelBlue, mark=none"))
	push!(p, Plots.Linear([-0.5, a, a], [fa, fa, -3], style="draw=none, mark=none, name path=A, forget plot"))
	push!(p, Plots.Linear([-0.5, b, b], [fb, fb, -3], style="draw=none, mark=none, name path=B, forget plot"))
	push!(p, Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];"))
	_μ = (fa+fb)/2
	σ = optimize(σ->(quantile(Normal(_μ, σ), 0.5+δ) - fb)^2, 0.1, 10.0).minimizer
	Q = Normal(_μ, σ)
	y_arr = range(fa, stop=fb, length=101)
	x_arr = [3*(pdf(Q, y)-pdf(Q,fa))-0.5 for y in y_arr]
	push!(p, Plots.Linear(x_arr, y_arr, style="solid, pastelBlue, mark=none"))
	push!(p, Plots.Node("sensitive",0,fb, style="above right"))

	Axis(p, width="9cm", xmin=-0.5, ymin=ymin, xlabel=L"x", ylabel=L"y", style="axis lines=left, axis on top, xtick=\\empty, ytick=\\empty")
end
plot(p)
# -

p = let
	_μ = x -> x^2 + 4/(1+abs(x))
	_σ = x -> sqrt(8/(1+abs(x)))
	Axis([
	    Plots.Linear(x -> _μ(x) + _σ(x), (-3,3), xbins=51, style="draw=none, mark=none, name path=A, forget plot"),
	    Plots.Linear(x -> _μ(x) - _σ(x), (-3,3), xbins=51, style="draw=none, mark=none, name path=B, forget plot"),
	    Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{\$\\pm \\sigma\$}"),
	    Plots.Linear(x -> _μ(x), (-3,3), xbins=51, style="solid, black, mark=none", legendentry="expected value"),
	], xmin=-3, xmax=3, width="10cm", height="5cm", xlabel=L"x", ylabel=L"y", style="legend style={draw=none, at={(0.5,-0.25)}, anchor=north, legend columns=-1, /tikz/every even column/.append style={column sep=0.5cm}}")
end
plot(p)

p = let
	_μ = x -> x^2 + 4/(1+abs(x))
	_σ = x -> sqrt(8/(1+abs(x)))
	p = Plots.Plot[]
	p_scatter = Plots.Plot[]
	for α in range(0.0,stop=1.0,length=11)
	    f = x -> α*_μ(x) + (1-α)*_σ(x)
	    xopt = 0.0
	    for x in range(-3,stop=0,length=101)
	        if f(x) < f(xopt)
	            xopt = x
	        end
	    end
	    push!(p, Plots.Linear(f, (-3,3), xbins=151, style="solid, color=pastelBlue!$(round(Int, 80α+20))"))
	    push!(p_scatter, Plots.Scatter([xopt, -xopt], f.([xopt, -xopt]), style="solid, mark=*, mark size=0.75, mark options={draw=black, fill=black}"))
	end
	append!(p, p_scatter)
	Axis(p, xmin=-3, xmax=3, width="10cm", height="5cm", xlabel=L"x", ylabel=L"\alpha \mu + (1-\alpha)\sigma", style="axis on top, clip marker paths=true, legend pos = outer north east")
end
plot(p)

using Distributions
p = let
	xmax = 7.0
	ymax = 0.4

	G = Truncated(Gamma(2.0,1.0), 0.0, xmax)

	μ_ = quadgk(x -> x*pdf(G, x), 0.0, xmax)[1]

	α = 0.1
	value_at_risk(α) = minimum(1-cdf(G, t) <= α ? t : Inf for t in range(0.0, stop=7.0, length=201))
	VaR = value_at_risk(α)
	cond_value_at_risk(α) = (1/α)*quadgk(t -> value_at_risk(1-t), 1-α, 1-0.001)[1]
	CVaR = cond_value_at_risk(α)

	xq = quantile(G, 1-α)

	Axis([
	    Plots.Linear(z->pdf(G, z), (0, xmax), style="solid, pastelBlue, mark=none"),
	    Plots.Linear(z->pdf(G, z), (xq, xmax), style="solid, pastelBlue, ultra thick, mark=none"),
	    Plots.Linear([μ_,μ_], [0, ymax], style="solid, gray, mark=none"),
	    Plots.Linear([VaR,VaR], [0, ymax], style="solid, gray, mark=none"),
	    Plots.Linear([CVaR,CVaR], [0, ymax], style="solid, gray, mark=none"),
	    Plots.Node("top \$1-\\alpha\$ quantile", (xq+xmax)/2+0.35, 0.08),
	], xlabel=L"y", ylabel=L"p(y)", style="xtick={$μ_,$VaR, $CVaR}, xticklabels={expected value, VaR, CVaR}", ymin=0, ymax=ymax, xmin=0, xmax=xmax, width="12cm", height="4cm"
	)
end
plot(p)

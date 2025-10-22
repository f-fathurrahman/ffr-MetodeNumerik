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

# # Surrogate Optimization
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

include("gp.jl")

# +

p = let
	N = 4
	y_dom = (-2,2)
	G = GroupPlot(N, 1, groupStyle="horizontal sep=0.25cm, ylabels at=edge left",
	                    style="width=4.5cm, xlabel=\$x\$, ylabel=\$y\$, xtick=\\empty, ytick=\\empty, xmin=0, xmax=8, ymin=$(y_dom[1]), ymax=$(y_dom[2])")

	f_true = x -> sin(x) + 0.25*sin(x/2) - 0.3*sin(2x)
	x_train = [3.5, 5.5]
	x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	GP = GaussianProcess()
	for x in x_train
		push!(GP, [x], f_true(x))
	end

	for n in 1 : N
	    p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	    # predict
	    μₚ, νₚ = predict(GP, X)
	    σₚ = sqrt.(νₚ)
	    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	    upperConfidence = μₚ + 1.96*σₚ
	    lowerConfidence = μₚ - 1.96*σₚ

	    # prediction-based optimization
	    pred_μ = x -> predict(GP, Vector{Float64}[Float64[x]])[1][1]
	    x_slice = myopt(pred_μ)

	    # slice
	    y_dom = (-2,2)
	    p_sample = Plots.Scatter([x_slice], [f_true(x_slice)], style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="next sample")
	    Plots.Linear([x_slice, x_slice], [y_dom[1], y_dom[2]], style="pastelBlue!40, mark=none")

	    ax = Axis([p_true, p_pred_μ, p_train,  p_sample])
	    if n == N
	        ax.legendPos="outer north east"
	    else
	        push!(ax.plots, Plots.Command("\\legend{}"))
	    end
	    push!(G, ax)

	    push!(GP, [x_slice], f_true(x_slice))
	    push!(x_train, x_slice)
	    x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	    X = [[x] for x in collect(x_arr)]
	end

	G
end
plot(p)
# -

using Optim
p = let
	N = 4
	M = 2
	y_dom = (-2,2)
	G = GroupPlot(N, M, groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, ylabels at=edge left, xlabels at=edge bottom",
	                    style="axis on top=true, clip marker paths=true, width=4.5cm, xlabel=\$x\$, ylabel=\$y\$, xtick=\\empty, ytick=\\empty, xmin=0, xmax=8, ymin=$(y_dom[1]), ymax=$(y_dom[2])")

	f_true = x -> sin(x) + 0.25*sin(x/2) - 0.3*sin(2x)
	x_train = [3.5, 5.5]
	x_arr = unique(sort!(append!(collect(range(0, stop=8, length=101)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	GP = GaussianProcess()
	for x in x_train
		push!(GP, [x], f_true(x))
	end

	for n in 1 : N*M
	    p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="only marks, solid, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	    # predict
	    μₚ, νₚ = predict(GP, X)
	    σₚ = sqrt.(νₚ)
	    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	    upperConfidence = μₚ + 1.96*σₚ
	    lowerConfidence = μₚ - 1.96*σₚ
	    p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
		p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
		p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

	    # prediction-based optimization
	    pred_ν = x->predict(GP, Vector{Float64}[Float64[x]])[2][1]
	    x_slice = myopt(x->-pred_ν(x))

	    # slice
	    y_dom = (-2,2)
	    p_sample = Plots.Scatter([x_slice], [f_true(x_slice)], style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="sampled")
	    Plots.Linear([x_slice, x_slice], [y_dom[1], y_dom[2]], style="pastelBlue!40, mark=none")

	    ax = Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_pred_μ, p_train, p_sample])
	    if n == N
	        ax.legendPos="outer north east"
	    else
	        push!(ax.plots, Plots.Command("\\legend{}"))
	    end
	    push!(G, ax)

	    push!(GP, [x_slice], f_true(x_slice))
	    push!(x_train, x_slice)

	    i = findfirst(x_arr .> x_slice)
	    if min(abs(x_arr[i] - x_slice), abs(x_arr[i-1] - x_slice)) > 1e-4
	    	insert!(x_arr, i, x_slice)
	    end
	    X = [[x] for x in collect(x_arr)]
	end

	G
end
plot(p)

using Optim
p = let
	N = 4
	M = 2
	y_dom = (-2,2)
	G = GroupPlot(N, M, groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, ylabels at=edge left, xlabels at=edge bottom",
	                    style="axis on top=true, clip marker paths=true, width=4.5cm, xlabel=\$x\$, ylabel=\$y\$, xtick=\\empty, ytick=\\empty, xmin=0, xmax=8, ymin=$(y_dom[1]), ymax=$(y_dom[2])")


	f_true = x -> sin(x) + 0.25*sin(x/2) - 0.3*sin(2x)
	x_train = [3.5, 5.5]
	x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	GP = GaussianProcess()
	for x in x_train
		push!(GP, [x], f_true(x))
	end

	for n in 1 : N*M
	    p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	    # predict
	    μₚ, νₚ = predict(GP, X)
	    σₚ = sqrt.(νₚ)
	    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	    upperConfidence = μₚ + 1.96*σₚ
	    lowerConfidence = μₚ - 1.96*σₚ
	    p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
		p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
		p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

	    # prediction-based optimization
	    pred_μ = x->predict(GP, Vector{Float64}[Float64[x]])[1][1]
	    pred_ν = x->predict(GP, Vector{Float64}[Float64[x]])[2][1]
	    f_opt = x->pred_μ(x) - sqrt(pred_ν(x))
	    x_slice = myopt(f_opt)
	    p_opt = Plots.Linear(x_arr, f_opt.(x_arr), style="solid, thick, pastelRed, mark=none", legendentry="lower bound")
	    p_opt_vert = Plots.Linear([x_slice, x_slice], [f_opt(x_slice), f_true(x_slice)], style="solid, pastelRed!40, mark=none, forget plot")

	    # slice
	    y_dom = (-2,2)
	    p_sample = Plots.Scatter([x_slice], [f_true(x_slice)], style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="sampled")
	    Plots.Linear([x_slice, x_slice], [y_dom[1], y_dom[2]], style="pastelBlue!40, mark=none")

	    ax = Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_opt_vert, p_opt, p_pred_μ, p_train, p_sample])
	    if n == N
	        ax.legendPos="outer north east"
	    else
	        push!(ax.plots, Plots.Command("\\legend{}"))
	    end
	    push!(G, ax)

	    push!(GP, [x_slice], f_true(x_slice))
	    push!(x_train, x_slice)
	    x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	    X = [[x] for x in collect(x_arr)]
	end

	G
end
plot(p)

# +
using Random
using Vec

p = let

	GP = GaussianProcess()

	# create random "true" function
	Random.seed!(3)
	x_arr = range(0, stop=8, length=101)
	X = [[x] for x in collect(x_arr)]
	y = rand(GP, X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	for i in [5, 20, 35, 50]
		push!(GP, X[i], y[i])
	end
	p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	# predict
	μₚ, νₚ = predict(GP, X)
	σₚ = sqrt.(νₚ)
	p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	upperConfidence = μₚ + 1.96*σₚ
	lowerConfidence = μₚ - 1.96*σₚ

	# slice
	slice_index = 70
	y_dom = (-3,3)
	x_slice = X[slice_index][1]
	μ_slice = μₚ[slice_index]
	σ_slice = σₚ[slice_index]
	N = Normal(μ_slice, σ_slice)
	y_arr = collect(range(y_dom[1], stop=y_dom[2], length=101))
	x_arr = [x_slice + pdf(N, y) for y in y_arr]
	p_slice = Plots.Linear(x_arr, y_arr, style="pastelBlue!40, mark=none")
	p_slice_vert = Plots.Linear([x_slice, x_slice], [y_arr[1], y_arr[end]], style="black, mark=none")

	y_min = minimum(GP.y)
	x_min = GP.X[argmin(GP.y)][1]
	p_fill_lo = Plots.Linear([x_slice, x_slice], [y_arr[1], y_min], style="draw=none, mark=none, name path=C, forget plot")
	y_arr_hi = collect(range(y_dom[1], stop=y_min, length=51))
	x_arr_hi = [x_slice + pdf(N, y) for y in y_arr_hi]
	p_fill_hi = Plots.Linear(x_arr_hi, y_arr_hi, style="draw=none, mark=none, name path=D, forget plot")
	p_fill = Plots.Command("\\addplot[pastelBlue!40, forget plot] fill between[of=C and D];")
	p_flat = Plots.Linear([x_min, x_slice], [y_min, y_min], style="solid, pastelBlue!40, mark=none, forget plot")

	arow_best_so_far = Plots.Command("\\draw [->, thick, shorten >=5pt] (axis description cs:0.3,0.75) -- (axis cs:$(x_min),$(y_min));")
	node_best_so_far = Plots.Node("\\small best so far", 0.3, 0.75, style="anchor=south east", axis="axis description cs")

	A = VecE2(x_slice-1, y_min-1)
	arow_prob_of_imp = Plots.Command("\\draw [->, thick] (axis cs:$(A.x),$(A.y)) -- (axis cs:$(x_slice+0.08),$(A.y+0.3));")
	node_prob_of_imp = Plots.Node("\\small probability of improvement", A.x, A.y, style="anchor=east", axis="axis cs")

	B = A + VecE2(0, -0.75)
	arow_query_point = Plots.Command("\\draw [->, thick, shorten >=5pt] (axis cs:$(B.x),$(B.y)) -- (axis cs:$(x_slice),$(B.y));")
	node_query_point = Plots.Node("\\small query point", B.x, B.y, style="anchor=east", axis="axis cs")

	Axis([p_true, p_pred_μ, p_train, p_slice, p_slice_vert, p_flat, p_fill_lo, p_fill_hi, p_fill,
	    node_best_so_far, arow_best_so_far, node_prob_of_imp, arow_prob_of_imp, node_query_point, arow_query_point],
	    xmin=0, xmax=8, style="ytick=\\empty, xtick=\\empty",
	    ymin=y_dom[1], ymax=y_dom[2],
	    xlabel=L"x", ylabel=L"y", width="9cm", legendPos="outer north east",
	)
end

plot(p)

# +
prob_of_improvement(y_min, μ, ν) = isapprox(ν, 0, atol=1e-4) ? 0.0 : cdf(Normal(μ, sqrt(ν)), y_min)

p = let

	N = 4
	y_dom = (-2,2)
	G = GroupPlot(N, 2, groupStyle="horizontal sep=0.25cm, ylabels at=edge left, vertical sep=0.25cm, xlabels at=edge bottom",
	                    style="axis on top=true, clip marker paths=true, width=4.5cm, xlabel=\$x\$, xtick=\\empty, ytick=\\empty, xmin=0, xmax=8, ymin=$(y_dom[1]), ymax=$(y_dom[2])")

	f_true = x -> sin(x) + 0.25*sin(x/2) - 0.3*sin(2x)
	x_train = [3.5, 5.5]
	x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	GP = GaussianProcess()
	for x in x_train
		push!(GP, [x], f_true(x[1]))
	end

	axes_lo = Axis[]

	for n in 1 : N
	    p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	    # predict
	    μₚ, νₚ = predict(GP, X)
	    σₚ = sqrt.(νₚ)
	    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	    upperConfidence = μₚ + 1.96*σₚ
	    lowerConfidence = μₚ - 1.96*σₚ
	    p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
		p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
		p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

	    # prediction-based optimization
	    pred_μ = x->predict(GP, Vector{Float64}[Float64[x]])[1][1]
	    pred_ν = x->predict(GP, Vector{Float64}[Float64[x]])[2][1]
	    y_min = minimum(GP.y)
	    f_opt = x->prob_of_improvement(y_min, pred_μ(x), pred_ν(x))
	    x_slice = myopt(x->-f_opt(x))
	    p_opt = Plots.Linear(x_arr, f_opt.(x_arr), style="solid, pastelRed, mark=none")
	    p_opt_vert = Plots.Linear([x_slice, x_slice], [f_opt(x_slice), f_true(x_slice)], style="solid, pastelRed!40, mark=none")

	    # slice
	    y_dom = (-2,2)
	    p_sample = Plots.Scatter([x_slice], [f_true(x_slice)], style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="sampled")
	    Plots.Linear([x_slice, x_slice], [y_dom[1], y_dom[2]], style="pastelBlue!40, mark=none")

	    push!(axes_lo, Axis([p_opt,
	                Plots.Scatter([x_slice], [f_opt(x_slice)], style="axis on top, only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}")
	                ], ymin=0, ymax=1))

	    ax = Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_pred_μ, p_train, p_sample])
	    if n == N
	        ax.legendPos="outer north east"
	    else
	        push!(ax.plots, Plots.Command("\\legend{}"))
	    end
	    if n == 1
	        ax.ylabel=L"y"
	        axes_lo[1].ylabel=L"P[f(\vect x) < y_\text{min}]"
	        axes_lo[1].style = "ytick={0,1}"
	    end
	    push!(G, ax)
	    push!(GP, [x_slice], f_true(x_slice))
	    push!(x_train, x_slice)
	    x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	    X = [[x] for x in collect(x_arr)]
	end
	for ax in axes_lo
	    push!(G, ax)
	end

	G
end
plot(p)

# +

function expected_improvement(y_min, μ, ν)
	if isapprox(ν, 0.0, atol=1e-4)
		return 0.0
	end
	σ = sqrt(ν)
    tmp = (y_min - μ)/σ
    p_imp = prob_of_improvement(y_min, μ, ν)
    p_ymin = pdf(Normal(μ, σ), y_min)
    return (y_min - μ)*p_imp + σ^2*p_ymin
end

p = let

	N = 4
	y_dom = (-2,2)
	G = GroupPlot(N, 2, groupStyle="horizontal sep=0.25cm, ylabels at=edge left, vertical sep=0.25cm, xlabels at=edge bottom",
	                    style="axis on top, clip marker paths=true, width=4.5cm, xlabel=\$x\$, xtick=\\empty, ytick=\\empty, xmin=0, xmax=8, ymin=$(y_dom[1]), ymax=$(y_dom[2])")

	f_true = x -> sin(x) + 0.25*sin(x/2) - 0.3*sin(2x)
	x_train = [3.5, 5.5]
	x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	GP = GaussianProcess()
	for x in x_train
		push!(GP, [x], f_true(x))
	end

	axes_lo = Axis[]

	for n in 1 : N
	    p_train = Plots.Scatter([x[1] for x in GP.X], GP.y, style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}", legendentry="fit points")

	    # predict
	    μₚ, νₚ = predict(GP, X)
	    σₚ = sqrt.(νₚ)
	    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	    upperConfidence = μₚ + 1.96*σₚ
	    lowerConfidence = μₚ - 1.96*σₚ
	    p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
		p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
		p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

	    # prediction-based optimization
	    pred_μ = x->predict(GP, Vector{Float64}[Float64[x]])[1][1]
	    pred_ν = x->predict(GP, Vector{Float64}[Float64[x]])[2][1]
	    y_min = minimum(GP.y)
	    f_opt = x->expected_improvement(y_min, pred_μ(x), pred_ν(x) + 1e-8)
	    x_slice = myopt(x->-f_opt(x))
	    p_opt = Plots.Linear(x_arr, f_opt.(x_arr), style="solid, pastelRed, mark=none")
	    p_opt_vert = Plots.Linear([x_slice, x_slice], [f_opt(x_slice), f_true(x_slice)], style="solid, pastelRed!40, mark=none")

	    # slice
	    y_dom = (-2,2)
	    p_sample = Plots.Scatter([x_slice], [f_true(x_slice)], style="only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}", legendentry="sampled")
	    Plots.Linear([x_slice, x_slice], [y_dom[1], y_dom[2]], style="pastelBlue!40, mark=none")

	    push!(axes_lo, Axis([p_opt,
	                Plots.Scatter([x_slice], [f_opt(x_slice)], style="axis on top, only marks, mark=*, mark size=1, mark options={draw=pastelRed, fill=pastelRed}")
	                ], ymin=0, ymax=1))

	    ax = Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_pred_μ, p_train, p_sample])
	    if n == N
	        ax.legendPos="outer north east"
	    else
	        push!(ax.plots, Plots.Command("\\legend{}"))
	    end
	    if n == 1
	        ax.ylabel=L"y"
	        axes_lo[1].ylabel = L"E[I(y)]"
	        axes_lo[1].style = "ytick={0,1}"
	    end
	    push!(G, ax)
	    push!(GP, [x_slice], f_true(x_slice))
	    push!(x_train, x_slice)
	    x_arr = unique(sort!(append!(collect(range(0, stop=8, length=201)), x_train)))
	    X = [[x] for x in collect(x_arr)]
	end
	for ax in axes_lo
	    push!(G, ax)
	end

	G
end
plot(p)

# +

if !@isdefined(GaussianProcess)
	include("../chapter/gp.jl")
end

using Distributions
p = let
	Z = Normal(0.0, 1.0)
	dom=(-2,6)

	function f_(x::Float64, z::Float64=0.0)
	    A = Normal(1.0, 1.0)
	    B = Normal(3.0, 0.5)
	    C = Normal(5.0, 0.25)
	    return 1 - (2pdf(A, x+z) + pdf(B, x+z) + pdf(C, x+z))
	end

	y_max = 0.8
	lo, hi = -0.8, 1.0

	p = Plots.Plot[]

	# append!(p, plot_transparent_interval((-0.6635, 3.8425), (lo,hi), style="green!50"))
	# append!(p, plot_transparent_interval((4.4838, 5.5095), (lo,hi), style="green!50"))
	append!(p, plot_transparent_intervals([(-0.6635, 3.8425), (4.4838, 5.5095)], (lo,hi), "green", 0.5))

	push!(p, Plots.Linear(f_, dom, xbins=201, style="solid, black, mark=none"))
	push!(p, Plots.Linear([dom[1], dom[2]], [y_max, y_max], style="solid, gray, mark=none"))
	push!(p, Plots.Scatter([5.0], [f_(5.0)], style="only marks, black, mark=*, mark size=1, mark options={draw=black, fill=black}"))
	push!(p, Plots.Node(L"x^*", 5.0, f_(5.0), style="right"))

	Axis(p, width="9cm", xlabel=L"x", ylabel=L"y", ymin=lo, ymax=hi, style="enlarge x limits=0, xtick=\\empty, ytick={$(y_max)}, yticklabel={\$y_\\text{max}\$}, axis on top",
	)

end
plot(p)

# +
using Random

p = let
	P_safe = 0.9
	σ_safe = optimize(x->(cdf(Normal(0.0,1.0), x) - P_safe)^2, -5, 5).minimizer
	β = σ_safe*σ_safe
	y_max = 0.8

	function f_(x::Float64, z::Float64=0.0)
	    A = Normal(1.0, 1.0)
	    B = Normal(3.0, 0.5)
	    C = Normal(5.0, 0.25)
	    return 1 - (2pdf(A, x+z) + pdf(B, x+z) + pdf(C, x+z))
	end
	f_true = x -> f_(x,0.0)
	dom = (-2,6)
	GP = GaussianProcess(ν = 0.01)

	x_train = [0.5,1.0,2.0,3.0,3.5]
	x_arr = unique(sort!(append!(collect(range(dom[1], stop=dom[2], length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)
	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="true")

	# create training set
	Random.seed!(0)
	for x in x_train
	    push!(GP, [x], f_true(x) + randn()*sqrt(GP.ν))
	end
	p_train = plot_GP_data(GP, legendentry="fit points")

	# predict
	(μₚ, νₚ) = predict(GP, X)
	σₚ = sqrt.(νₚ)
	p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	upperConfidence = μₚ + σ_safe*σₚ
	lowerConfidence = μₚ - σ_safe*σₚ
	p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
	p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
	p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")
	p_lowerbound = Plots.Linear(x_arr, lowerConfidence, style="solid, thick, pastelGreen, mark=none", legendentry=L"\ell(x)")
	p_upperbound = Plots.Linear(x_arr, upperConfidence, style="solid, thick, pastelRed, mark=none", legendentry=L"u(x)")

	Axis([p_pred_hi, p_pred_lo, p_pred, p_true, p_pred_μ, p_lowerbound, p_upperbound, p_train],
				style="legend pos=outer north east, axis on top=true, clip marker paths=true, width=9cm, height=4cm, xlabel=\$x\$, xtick=\\empty, ytick=\\empty, enlarge x limits=0, legend style={draw=none}")
end
plot(p)

# +
using Random

p = let
	lo, hi = -1.0, 1.2

	P_safe = 0.9
	σ_safe = optimize(x->(cdf(Normal(0.0,1.0), x) - P_safe)^2, -5, 5).minimizer
	β = σ_safe*σ_safe
	y_max = 0.8

	function f_(x::Float64, z::Float64=0.0)
	    A = Normal(1.0, 1.0)
	    B = Normal(3.0, 0.5)
	    C = Normal(5.0, 0.25)
	    return 1 - (2pdf(A, x+z) + pdf(B, x+z) + pdf(C, x+z))
	end
	f_true = x -> f_(x,0.0)
	dom = (-2,6)
	GP = GaussianProcess(ν = 0.01)

	x_train = [0.5,1.0,2.0,3.0,3.5]
	x_arr = unique(sort!(append!(collect(range(dom[1], stop=dom[2], length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry=L"\text{objective function } f")

	# create training set
	Random.seed!(0)
	for x in x_train
	    push!(GP, [x], f_true(x) + randn()*sqrt(GP.ν))
	end
	p_train = plot_GP_data(GP, legendentry="fit points")

	# predict
	(μₚ, νₚ) = predict(GP, X)
	σₚ = sqrt.(νₚ)
	p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	upperConfidence = μₚ + σ_safe*σₚ
	lowerConfidence = μₚ - σ_safe*σₚ
	p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
	p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
	p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence interval}")
	p_lowerbound = Plots.Linear(x_arr, lowerConfidence, style="solid, thick, pastelGreen, mark=none", legendentry=L"\ell(x)")
	p_upperbound = Plots.Linear(x_arr, upperConfidence, style="solid, thick, pastelRed, mark=none", legendentry=L"u(x)")

	p = Plots.Plot[]
	push!(p, p_pred_hi)
	push!(p, p_pred_lo)
	push!(p, p_pred)
	push!(p, p_true)
	push!(p, p_train)
	push!(p, Plots.Linear([dom[1], dom[2]], [y_max, y_max], style="solid, thick, gray, mark=none", legendentry="safety threshold"))

	safe_regions = get_safe_regions(GP, β, dom..., y_max)
	append!(p, plot_transparent_intervals(safe_regions, (lo,hi), "green", 0.5, "estimated safe region \$\\mathcal{S}\$"))

	G = GroupPlot(1, 2, groupStyle="vertical sep=0.5cm, xlabels at=edge bottom",
	    style="axis on top=true, width=9cm, height=4cm, xlabel=\$x\$, xtick=\\empty, enlarge x limits=0")
	push!(G, Axis(p, xlabel=L"x", ylabel=L"y", ymin=lo, ymax=hi, style="legend pos=outer north east, ytick=\\empty, legend style={draw=none}"))

	push!(G, Axis(
	    [
	    Plots.Linear(x->prob_is_safe(predict(GP, [x]), y_max), dom, xbins=300, style="solid, green, thick, mark=none"),
	    Plots.Linear([dom[1], dom[2]], [P_safe, P_safe], style="solid, gray, mark=none"),
	    Plots.Node(L"P_\text{safe}", 3, P_safe, style="below")
	    ], ymin=0, ymax=1, ylabel="safety probability", style="enlarge x limits=0, ytick={0,1} "
	))
end
plot(p)

# +
using Random

p = let
	lo, hi = -1.0, 1.2

	P_safe = 0.9
	σ_safe = optimize(x->(cdf(Normal(0.0,1.0), x) - P_safe)^2, -5, 5).minimizer
	β = σ_safe*σ_safe
	y_max = 0.8

	function f_(x::Float64, z::Float64=0.0)
	    A = Normal(1.0, 1.0)
	    B = Normal(3.0, 0.5)
	    C = Normal(5.0, 0.25)
	    return 1 - (2pdf(A, x+z) + pdf(B, x+z) + pdf(C, x+z))
	end
	f_true = x -> f_(x,0.0)
	dom = (-2,6)
	GP = GaussianProcess(ν = 0.01)

	x_train = [0.5,1.0,2.0,3.0,3.5]
	x_arr = unique(sort!(append!(collect(range(dom[1], stop=dom[2], length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="objective function")

	# create training set
	Random.seed!(0)
	for x in x_train
	    push!(GP, [x], f_true(x) + randn()*sqrt(GP.ν))
	end
	p_train = plot_GP_data(GP, legendentry="fit points")

	# predict
	(μₚ, νₚ) = predict(GP, X)
	σₚ = sqrt.(νₚ)
	p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
	upperConfidence = μₚ + σ_safe*σₚ
	lowerConfidence = μₚ - σ_safe*σₚ
	p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
	p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
	p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

	p = Plots.Plot[]
	push!(p, p_pred_hi)
	push!(p, p_pred_lo)
	push!(p, p_pred)
	push!(p, p_true)
	push!(p, p_train)
	push!(p, Plots.Linear([dom[1], dom[2]], [y_max, y_max], style="solid, thick, gray, mark=none", legendentry="safety threshold"))

	(best_hi, i) = findmin(upperConfidence)
	push!(p, Plots.Linear([dom[1], dom[2]], [best_hi, best_hi], style="solid, thick, pastelRed, mark=none", legendentry="best upper bound"))

	ax3 = Axis(p, xlabel=L"x", ylabel=L"y", ymin=lo, ymax=hi, style="legend pos=outer north east, ytick=\\empty, legend style={draw=none}, axis on top=true, width=9.5cm, height=4cm, xlabel=\$x\$, xtick=\\empty, enlarge x limits=0")

	safe_regions = get_safe_regions(GP, β, dom..., y_max)
	ax1 = Axis(plot_transparent_intervals(safe_regions, (lo,hi), "green", 0.5, "estimated safe region \$\\mathcal{S}\$"),
	        xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi, style="ytick=\\empty, xtick=\\empty, legend pos=outer north east, axis lines=none, legend style={draw=none}, width=9.5cm, height=2.05cm")

	M_regions = get_potential_maximizers(GP, β, safe_regions)
	ax2 = Axis(plot_transparent_intervals(M_regions, (lo,hi), "magenta", 0.5, "potential minimizers \$\\mathcal{M}\$"),
	        xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi, style="ytick=\\empty, xtick=\\empty, legend pos=outer north east, axis lines=none, legend style={draw=none}, width=9.5cm, height=2.05cm")

	G = GroupPlot(1, 3, groupStyle="vertical sep=0cm")
	push!(G, ax1)
	push!(G, ax2)
	push!(G, ax3)
	G
end
plot(p)

# +
using Random

p = let
    lo, hi = -1.0, 1.2

    P_safe = 0.9
    σ_safe = optimize(x->(cdf(Normal(0.0,1.0), x) - P_safe)^2, -5, 5).minimizer
    β = σ_safe*σ_safe
    y_max = 0.8

    function f_(x::Float64, z::Float64=0.0)
        A = Normal(1.0, 1.0)
        B = Normal(3.0, 0.5)
        C = Normal(5.0, 0.25)
        return 1 - (2pdf(A, x+z) + pdf(B, x+z) + pdf(C, x+z))
    end
    f_true = x -> f_(x,0.0)
    dom = (-2,6)
    GP = GaussianProcess(ν = 0.01)

    lo, hi = -1.0, 1.2

    x_train = [0.5,1.0,2.0,3.0,3.5]
    x_arr = unique(sort!(append!(collect(range(dom[1], stop=dom[2], length=201)), x_train)))
    X = [[x] for x in collect(x_arr)]
    y = (x -> f_true(x[1])).(X)

    p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none", legendentry="objective function")

    # create training set
    Random.seed!(0)
    GP = GaussianProcess(ν = 0.01)
    for x in x_train
        push!(GP, [x], f_true(x) + randn()*sqrt(GP.ν))
    end
    p_train = plot_GP_data(GP, legendentry="fit points")

    # predict
    (μₚ, νₚ) = predict(GP, X)
    σₚ = sqrt.(νₚ)
    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none", legendentry="predicted")
    upperConfidence = μₚ + σ_safe*σₚ
    lowerConfidence = μₚ - σ_safe*σₚ
    p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
    p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
    p_pred = Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B];\\addlegendentry{confidence region}")

    p = Plots.Plot[]
    push!(p, p_pred_hi)
    push!(p, p_pred_lo)
    push!(p, p_pred)
    push!(p, p_true)
    push!(p, p_train)
    push!(p, Plots.Linear([dom[1], dom[2]], [y_max, y_max], style="solid, thick, gray, mark=none", legendentry="safety threshold"))

    ax3 = Axis(p, xlabel=L"x", ylabel=L"y", ymin=lo, ymax=hi, style="legend pos=outer north east, ytick=\\empty, legend style={draw=none}, axis on top=true, width=9.5cm, height=4cm, xlabel=\$x\$, xtick=\\empty, enlarge x limits=0")

    safe_regions = get_safe_regions(GP, β, dom..., y_max)
    ax1 = Axis(plot_transparent_intervals(safe_regions, (lo,hi), "green", 0.5, "estimated safe region \$\\mathcal{S}\$"),
            xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi, style="ytick=\\empty, xtick=\\empty, legend pos=outer north east, axis lines=none, legend style={draw=none}, width=9.5cm, height=2.05cm")

    d = (x,x′) -> norm(x-x′,2)
    L = 4.0

    E_regions = get_potential_expanders(GP, β, safe_regions, L, d, y_max)
    ax2 = Axis(plot_transparent_intervals(E_regions, (lo,hi), "orange", 0.5, "expected expanders \$\\mathcal{E}\$"),
            xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi, style="ytick=\\empty, xtick=\\empty, legend pos=outer north east, axis lines=none, legend style={draw=none}, width=9.5cm, height=2.05cm")

    G = GroupPlot(1, 3, groupStyle="vertical sep=0cm")
    push!(G, ax1)
    push!(G, ax2)
    push!(G, ax3)
    G
end
plot(p)

# +
using Random

p = let
	M = 9
	G = GroupPlot(1, 4M, groupStyle="vertical sep=0cm, xlabels at=edge bottom, xticklabels at=edge bottom", style="xlabel=\$x\$, width=10cm")
	region_height = "1.75cm"

	P_safe = 0.9
	σ_safe = optimize(x->(cdf(Normal(0.0,1.0), x) - P_safe)^2, -5, 5).minimizer
	β = σ_safe*σ_safe
	y_max = 0.8

	function f_(x::Float64, z::Float64=0.0)
	    A = Normal(1.0, 1.0)
	    B = Normal(3.0, 0.5)
	    C = Normal(5.0, 0.25)
	    return 1 - (2pdf(A, x+z) + pdf(B, x+z) + pdf(C, x+z))
	end
	f_true = x -> f_(x,0.0)
	dom = (-2,6)
	GP = GaussianProcess(ν = 0.01)

	lo, hi = -1.0, 1.2

	d = (x,x′) -> norm(x-x′,2)
	L = 4.0

	x_train = [2.5]
	x_arr = unique(sort!(append!(collect(range(dom[1], stop=dom[2], length=201)), x_train)))
	X = [[x] for x in collect(x_arr)]
	y = (x -> f_true(x[1])).(X)

	p_true = Plots.Linear(x_arr, y, style="solid, black, mark=none")

	# create training set
	Random.seed!(0)
	for x in x_train
	    push!(GP, [x], f_true(x) + randn()*sqrt(GP.ν))
	end

	for i in 1 : M
	    p_train = plot_GP_data(GP)

	    # predict
	    (μₚ, νₚ) = predict(GP, X)
	    σₚ = sqrt.(νₚ)
	    p_pred_μ = Plots.Linear(x_arr, μₚ, style="solid, thick, pastelBlue, mark=none")
	    upperConfidence = μₚ + σ_safe*σₚ
	    lowerConfidence = μₚ - σ_safe*σₚ
	    p_pred_hi = Plots.Linear(x_arr, upperConfidence, style="draw=none, mark=none, name path=A, forget plot")
	    p_pred_lo = Plots.Linear(x_arr, lowerConfidence, style="draw=none, mark=none, name path=B, forget plot")
	    p_pred = Plots.Command("\\addplot[pastelBlue!40, forget plot] fill between[of=A and B];")

	    p = Plots.Plot[]
	    push!(p, p_pred_hi)
	    push!(p, p_pred_lo)
	    push!(p, p_pred)
	    push!(p, p_true)
	    push!(p, p_train)
	    push!(p, Plots.Linear([dom[1], dom[2]], [y_max, y_max], style="solid, gray, mark=none"))

	    if i == 1
	        safe_regions = get_safe_regions(GP, β, dom..., y_max)
	        push!(G, Axis(plot_transparent_intervals(safe_regions, (lo,hi), "green", 0.5, "estimated safe region"),
	                      xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi,
	                      style="ytick=\\empty, xtick=\\empty, axis lines=none, height=$region_height, legend style={draw=none, at={(1.01,1)}, anchor=north west}"))

	        M_regions = get_potential_maximizers(GP, β, safe_regions)
	        push!(G, Axis(plot_transparent_intervals(M_regions, (lo,hi), "magenta", 0.5, "potential minimizers"),
	                      xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi,
	                      style="ytick=\\empty, xtick=\\empty, axis lines=none, height=$region_height, legend style={draw=none, at={(1.01,0.5)}, anchor=north west}"))

	        E_regions = get_potential_expanders(GP, β, safe_regions, L, d, y_max)
	        push!(G, Axis(plot_transparent_intervals(E_regions, (lo,hi), "orange", 0.5, "potential expanders"),
	                      xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi,
	                      style="ytick=\\empty, xtick=\\empty, axis lines=none, height=$region_height, legend style={draw=none, at={(1.01,0)}, anchor=north west}"))
	    else
	        safe_regions = get_safe_regions(GP, β, dom..., y_max)
	        push!(G, Axis(plot_transparent_intervals(safe_regions, (lo,hi), "green", 0.5),
	                      xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi,
	                      style="ytick=\\empty, xtick=\\empty, axis lines=none, height=$region_height"))

	        M_regions = get_potential_maximizers(GP, β, safe_regions)
	        push!(G, Axis(plot_transparent_intervals(M_regions, (lo,hi), "magenta", 0.5),
	                      xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi,
	                      style="ytick=\\empty, xtick=\\empty, axis lines=none, height=$region_height"))

	        E_regions = get_potential_expanders(GP, β, safe_regions, L, d, y_max)
	        push!(G, Axis(plot_transparent_intervals(E_regions, (lo,hi), "orange", 0.5),
	                      xmin=dom[1], xmax=dom[2], ymin=lo, ymax=hi,
	                      style="ytick=\\empty, xtick=\\empty, axis lines=none, height=$region_height"))
	    end

	    x_sample = [NaN]
	    w_best = -Inf
	    Δx = 0.01
	    for (a,b) in vcat(M_regions, E_regions)
	        m = ceil(Int, (b-a)/Δx)
	        for x in range(a, stop=b, length=m)
	            w = width(GP, [x], β)
	            if w > w_best
	                w_best = w
	                x_sample[1] = x
	            end
	        end
	    end
	    push!(p, Plots.Scatter(x_sample, [f_true(x_sample[1])], style="only marks, black, mark=o, mark size=1.5, mark options={draw=black}"))

	    push!(GP, x_sample, f_true(x_sample[1]) + randn()*sqrt(GP.ν))

	    ax = Axis(p, ylabel=L"y", ymin=lo, ymax=hi, style="ytick=\\empty, axis on top=true, height=2.5cm, xtick=\\empty, enlarge x limits=0")
	    push!(G, ax)
	end
	G
end
plot(p)

# +
using Colors
using Random

p = let
	G = GroupPlot(4, 4, groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, ylabels at=edge left, xlabels at=edge bottom, xticklabels at=edge bottom, yticklabels at=edge left",
	    style="xlabel=\$x_1\$, ylabel=\$x_2\$, width=4cm, height=4cm, contour/labels=false, view={0}{90}"
	)

	function flower(x; a=1, b=1, c=4)::Float64
	    if isapprox(norm(x), 0.0)
	        return 0.0
	    end
		return a*norm(x) + b*sin(c*atan(x[2], x[1]))
	end
	f_(x::Vector{Float64}) = flower(x)
	f_(x1::Float64, x2::Float64) = f_([x1,x2])

	xdomain = ( -3, 3)
	ydomain = ( -3, 3)

	y_max = 2.0

	function plot_GP_mean(GP, β)
	    m, k, ν = GP.m, GP.k, GP.ν
	    iK = inv(K(GP.X, GP.X, k) + ν*I)

	    u_pred = (x,y) -> begin
	                X_pred = Vector{Float64}[[x,y]]
	                tmp = K(X_pred, GP.X, k) * iK
	                μₚ = μ(X_pred, m) + tmp*(GP.y - μ(GP.X, m))
	                νₚ = diag(K(X_pred, X_pred, k) + ν*I - tmp*K(GP.X, X_pred, k))
	                return μₚ[1] + sqrt(β*νₚ[1])
	            end

	    plots = Plots.Plot[]
	    push!(plots, Plots.Image(u_pred, xdomain, ydomain, xbins=600, ybins=600, colormap=pasteljet, colorbar = false))
	    push!(plots, Plots.Contour(u_pred, xdomain, ydomain, levels=[y_max], contour_style="draw color=white", style="white", xbins=101, ybins=101))
	    push!(plots, Plots.Scatter([x[1] for x in GP.X], [x[2] for x in GP.X], style="only marks, white, mark=*, mark size=0.75, mark options={draw=white, fill=white}"))

	    Axis(plots, xmin=xdomain[1], xmax=xdomain[2], ymin=ydomain[1], ymax=ydomain[2])
	end

	GP = GaussianProcess(m=x->y_max + 0.5, ν = 0.01)
	β = 10.0
	n = 51
	X = Array{Vector{Float64}}(undef, n*n)
	i = 0
	for x1 in range(xdomain[1], stop=xdomain[2], length=n)
	    for x2 in range(ydomain[1], stop=ydomain[2], length=n)
	        X[i+=1] = [x1,x2]
	    end
	end
	i = argmin([norm([-2,1]-x,2) for x in X])

	m = length(X)
	u, ℓ = fill(Inf, m), fill(-Inf, m)
	S, M, E = falses(m), falses(m), falses(m)

	Random.seed!(0)
	push!(GP, X[i], f_(X[i]) + randn()*GP.ν)
	update_confidence_intervals!(GP, X, u, ℓ, β);
	S[:] = u .≤ y_max
	best_val, i_best = findmin(u[S])
	i_best = something(findfirst(isequal(i_best), cumsum(S)), 0)

	ax = plot_GP_mean(GP, β)
	push!(ax, Plots.Scatter([X[i_best][1]], [X[i_best][2]], style="only marks, white, mark=+, mark size=3, mark options={draw=white, fill=white}"))
	push!(G, ax)

	for k in 2 : prod(G.dimensions)
	    compute_sets!(GP, S, M, E, X, u, ℓ, y_max, β)
	    i = get_new_query_point(M, E, u, ℓ)

	    push!(GP, X[i], f_(X[i]))
	    update_confidence_intervals!(GP, X, u, ℓ, β);
	    S[:] = u .≤ y_max
	    best_val, i_best = findmin(u[S])
	    i_best = something(findfirst(isequal(i_best), cumsum(S)), 0)


	    ax = plot_GP_mean(GP, β)
	    push!(ax, Plots.Scatter([X[i_best][1]], [X[i_best][2]], style="only marks, white, mark=+, mark size=3, mark options={draw=white, fill=white}"))
	    push!(G, ax)
	end

	G
end
plot(p)

# +
p = let
function flower(x; a=1, b=1, c=4)::Float64
	if isapprox(norm(x), 0.0)
		return 0.0
	end
	return a*norm(x) + b*sin(c*atan(x[2], x[1]))
end

y_max = 2.0
xdomain = ( -3, 3)
ydomain = ( -3, 3)

plots = Plots.Plot[]
push!(plots, Plots.Image((x,y)->flower([x,y]), xdomain, ydomain, xbins=600, ybins=600, colormap=pasteljet, colorbar = false))
push!(plots, Plots.Contour((x,y)->flower([x,y]), xdomain, ydomain, levels=[y_max], contour_style="draw color=white", style="white", xbins=101, ybins=101))
Axis(plots, width="4cm", height="4cm", xlabel=L"x_1", ylabel=L"x_2", style="contour/labels=false, axis equal, view={0}{90}")
end
plot(p)
# -



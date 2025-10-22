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

# # Direct
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Optim

p = let

	f = x -> x[1]^2 - 0.9*x[1]*x[2] + x[2]^2

	xdomain = (-2, 1)
	ydomain = (-2, 1)

	coordinate = 1

	x0 = [-1.5, -1.5]
	pts = Vector{Float64}[x0]
	dir = zeros(Float64, length(x0))
	for i in 1 : 10
		x = pts[end]
		dir[coordinate] = 1.0
		α = optimize(α->f(x + α*dir), -100.0, 100.0).minimizer
		push!(pts, x + α*dir)
		dir[coordinate] = 0.0
		coordinate = mod(coordinate,length(x))+1
	end

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[0.1,0.2,0.5,1,2,3,4], style="width=\\columnwidth"))
	push!(plots, Plots.Linear3([x[1] for x in pts], [x[2] for x in pts], [f(x) for x in pts], style="black, solid, thick, mark=none"))
	Axis(plots, width="1.1*9cm", height="1.1*9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)
# -


    p = let
        f = x -> abs(x[1]+x[2]) + abs(x[2] - x[1]) - 7exp(-(x[1]-1.5)^2 - (x[2]-1.5)^2)*(x[1]+x[2])

        xdomain = (-1, 3)
        ydomain = (-1, 3)

        plots = Plots.Plot[]
        push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[-15,-10,-5,-1,0,0.5,1,1.5,2,2.5,3], style="width=\\columnwidth"))
        push!(plots, Plots.Linear3([0,0,0.5], [0.5,0,0], [f([0,0.5]), f([0,0]), f([0.5,0])], style="black, solid, mark=none"))
        Axis(plots, width="1.1*9cm", height="1.1*9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
    end

    plot(p)

# +
using Optim

p = let

	basis(i, n) = [k == i ? 1. : 0. for k in 1 : n]

	f = x -> x[1]^2 - 0.9*x[1]*x[2] + x[2]^2

	xdomain = (-2, 1)
	ydomain = (-2, 1)

	x0 = [-1.5, -1.5]
	n = length(x0)
	pts_vanilla = Vector{Float64}[x0]
	pts_accel = deepcopy(pts_vanilla)
	for i in 1 : 6
		x = pts_vanilla[end]
		d = basis(mod1(i,n), n)
		α = optimize(α->f(x + α*d), -100.0, 100.0).minimizer
		push!(pts_vanilla, x + α*d)

		x = pts_accel[end]
		d = basis(mod1(i,n), n)
		α = optimize(α->f(x + α*d), -100.0, 100.0).minimizer
		push!(pts_accel, x + α*d)
		if mod1(i,n) == n
			x = pts_accel[end]
			d = x - pts_accel[end-n]
			α = optimize(α->f(x + α*d), -100.0, 100.0).minimizer
			push!(pts_accel, x + α*d)
		end
	end
	pts_accel = pts_accel[1:length(pts_vanilla)]

	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[0.1,0.2,0.5,1,2,3,4], style="width=\\columnwidth, forget plot"))
	push!(plots, Plots.Linear3([x[1] for x in pts_vanilla],
	                           [x[2] for x in pts_vanilla],
	                           [f(x) for x in pts_vanilla],
	                           style="pastelBlue, solid, thick, mark=none",
	                           legendentry="original"))
	push!(plots, Plots.Linear3([x[1] for x in pts_accel],
	                           [x[2] for x in pts_accel],
	                           [f(x) for x in pts_accel],
	                           style="pastelRed, solid, thick, mark=none",
	                           legendentry="accelerated"))
	Axis(plots, width="1.1*9cm", height="1.1*9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}, legend cell align=left, legend style={at={(0.5,-0.12)}, anchor=north, legend columns=1},")
end

plot(p)

# +
using Optim
p = let
	function basis(i, n)
	    eᵢ = zeros(n)
	    eᵢ[i] = 1
	    eᵢ
	end

	f = x -> x[1]^2 - 0.9*x[1]*x[2] + x[2]^2

	xdomain = (-2, 1)
	ydomain = (-2, 1)


	x0 = [-1.5, -1.5]
	pts = Vector{Float64}[x0]

	n = length(x0)
	U = [basis(i,n) for i in 1 : n]

	for i in 1 : 3
		x = pts[end]
		for i in 1 : n
	        dir = U[i]
	        α = optimize(α->f(x + α*dir), -100.0, 100.0).minimizer
	        x += α*dir
	        push!(pts, x)
	    end
	    for i in 1 : n-1
	        U[i] = U[i+1]
	    end
	    U[n] = dir = x - x0
	    α = optimize(α->f(x + α*dir), -100.0, 100.0).minimizer
	    x0 = x + α*dir
	    push!(pts, x0)
	end


	plots = Plots.Plot[]
	push!(plots, Plots.Contour(f, xdomain, ydomain, levels=[0.1,0.2,0.5,1,2,3,4], style="width=\\columnwidth"))
	push!(plots, Plots.Linear3([x[1] for x in pts], [x[2] for x in pts], [f(x) for x in pts], style="black, thick, solid, mark=none"))
	Axis(plots, width="1.1*9cm", height="1.1*9cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
end

plot(p)
# -


	p = let
		g = GroupPlot(4,1,groupStyle="xlabels at=edge bottom, ylabels at=edge left, horizontal sep=0.5cm")

		f = x -> -exp(-(x[1]*x[2] - 1.5)^2 -(x[2]-1.5)^2)

		basis(i, n) = [k == i ? 1. : 0. for k in 1 : n]

		xdomain = (0, 3)
		ydomain = (0, 3)

		p_contour = Plots.Contour(f, xdomain, ydomain, style="width=\\columnwidth", xbins=151, ybins=151) # , levels=[-15,-10,-5,-1,0,0.5,1,1.5,2,2.5,3]

		function add_hooke_jeeves_pts!(plots, x, α)
			push!(plots, Plots.Linear3([x[1] - α,x[1] + α,x[1],x[1]], [x[2],x[2],x[2]-α, x[2]+α], [f([x[1]-α,x[2]]), f([x[1]+α,x[2]]), f([x[1],x[2]-α]), f([x[1],x[2]+α])], style="only marks, mark=*, mark size=1, mark options={draw=black, fill=black}"))
			push!(plots, Plots.Linear3([x[1]], [x[2]], [f([x[1],x[2]])], style="only marks, mark=*, mark size=1.75, mark options={draw=black, fill=black}"))
		end

		x = [0.75, 0.75]
		α = 0.5
		γ = 0.5
		y, n = f(x), length(x)
		for i in 1 : 4

			plots = Plots.Plot[]
			push!(plots, p_contour)
			add_hooke_jeeves_pts!(plots, x, α)
			ax = Axis(plots, ymin=0, xmin=0, ymax=3, xmax=3, width="5cm", height="5cm", xlabel=L"x_1", ylabel=L"x_2", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
			push!(g, ax)

	        improved = false
	        x_best, y_best = x, y
	        for i in 1 : n
            	for sgn in (-1,1)
    	            x′ = x + sgn*α*basis(i, n)
    	            y′ = f(x′)
    	            if y′ < y_best
    	                x_best, y_best, improved = x′, y′, true
    	            end
                end
	        end
            x, y = x_best, y_best

	        if !improved
	            α *= γ
	        end
	    end
		g
	end

	plot(p)

		using Statistics

		p = let

			g = GroupPlot(4,3,groupStyle="xlabels at=edge bottom, ylabels at=edge left, horizontal sep=0.25cm, vertical sep=0.25cm")

			f = x -> -exp(-(x[1]*x[2] - 1.5)^2 -(x[2]-1.5)^2)

			xdomain = (0, 3)
			ydomain = (0, 3)

			p_contour = Plots.Contour(f, xdomain, ydomain, xbins=151, ybins=151)

			function add_simplex!(S_hist)

			    plots = Plots.Plot[]
			    push!(plots, p_contour)

			    for (i,S) in enumerate(S_hist)
			        Scirc = push!(deepcopy(S), S[1])
			        push!(plots, Plots.Linear3([x[1] for x in Scirc], [x[2] for x in Scirc], [f(x) for x in Scirc], style="solid, black, mark=none, line join=round, opacity=$(0.7^(length(S_hist) - i))")) # $
			    end
			    ax = Axis(plots, ymin=0, xmin=0, ymax=3, xmax=3, width="5.25cm", height="5.25cm", style="xtick=\\empty, ytick=\\empty, contour/labels=false, axis equal, view={0}{90}")
			    push!(g, ax)
			end

			α = 1.0
			β = 2.0
            γ = 0.5

			S = [[0.25,0.25], [0.75,0.5], [0.75,1.25]]
			S_hist = [deepcopy(S)]

			delta = Inf
			Y = [f(x) for x in S]

			add_simplex!(S_hist)

			for iter in 1 : g.dimensions[1]*g.dimensions[2]-1
			    # sort lowest to highest
			    p = sortperm(Y)
			    S, Y = S[p], Y[p]
			    xl, yl = S[1], Y[1] # lowest
			    xh, yh = S[end], Y[end] # highest
			    xs, ys = S[end-1], Y[end-1] # second-highest
			    xm = mean(S[1:end-1]) # centroid
			    xr = xm + α*(xm - xh)   # reflection point
			    yr = f(xr)

			    if yr < yl
			        xe = xm + β*(xr-xm)
			        ye = f(xe)
			        S[end], Y[end] = ye < yr ? (xe, ye) : (xr, yr)
			    elseif yr ≥ ys
			        if yr < yh
			            xh, yh, S[end], Y[end] = xr, yr, xr, yr
			        end
			        xc = xm + γ*(xh - xm)
			        yc = f(xc)
			        if yc > yh
			            for i in 2 : length(Y)
			                S[i] = (S[i] + xl)/2
			                Y[i] = f(S[i])
			            end
			        else
			            S[end], Y[end] = xc, yc
			        end
			    else
			        S[end], Y[end] = xr, yr
			    end

			    push!(S_hist, deepcopy(S))
			    add_simplex!(S_hist)
			end

			g
		end

		plot(p)

# +
using LinearAlgebra

struct FuncEval
    x::Vector{Float64}
    y::Float64
end

struct Cell2D
    x::Vector{Float64}
    y::Float64
    w::Vector{Float64}
end

p = let

	xdomain = (-2.0, 2.0)
	ydomain = (-2.0, 2.0)
	δ = (xdomain[2] - xdomain[1])/6

	fevals = [
	    FuncEval([ 0.0, 0.0], 0.0),
	    FuncEval([ 0.0, 2δ],  1.0),
	    FuncEval([ 0.0,-2δ],  0.0),
	    FuncEval([ 2δ,  0.0],-0.5),
	    FuncEval([-2δ,  0.0], 0.2),
	]

	function lowerbound(x::Vector{Float64}, L::Float64, fevals::Vector{FuncEval})
	    maximum(fe.y - L*norm(x - fe.x, 2) for fe in fevals)
	end

	G = GroupPlot(2,1,groupStyle="horizontal sep=1cm, ylabels at=edge left, yticklabels at=edge left",
				  style="width=6.5cm, height=6.5cm, xlabel=\$x_1\$, ylabel=\$x_2\$, axis equal, view={0}{90}, every tick/.style={white}")

	plots = Plots.Plot[]
	push!(plots, Plots.Image((a,b)->lowerbound([a,b], 1.0, fevals), xdomain, ydomain, xbins=600, ybins=600, colormap = pasteljet, colorbar = false))
	push!(plots, Plots.Scatter([fe.x[1] for fe in fevals], [fe.x[2] for fe in fevals], style="only marks, mark=*, mark size=1, mark options={draw=white, fill=white}"))
	push!(G, Axis(plots, title="Lipschitz lower bound"))

	cells = [
	    Cell2D([ 0.0, 0.0], 0.0, [4/3,4/3]),
	    Cell2D([ 0.0, 2δ],  0.5, [4,  4/3]),
	    Cell2D([ 0.0,-2δ],  0.0, [4,  4/3]),
	    Cell2D([ 2δ,  0.0],-0.5, [4/3,4/3]),
	    Cell2D([-2δ,  0.0], 0.2, [4/3,4/3]),
	]

	function cell_lower_bound(x, L, cells)
	    for cell in cells
	        Δ = abs.(x - cell.x)
	        if all(Δ .<= cell.w./2)
	            return cell.y - L*norm(Δ, 2)
	        end
	    end
	    return 0.0
	end

	plots = Plots.Plot[]
	push!(plots, Plots.Image((a,b)->cell_lower_bound([a,b], 1.0, cells), xdomain, ydomain, xbins=600, ybins=600, colormap = pasteljet, colorbar = false))
	push!(plots, Plots.Scatter([fe.x[1] for fe in fevals], [fe.x[2] for fe in fevals], style="only marks, mark=*, mark size=1, mark options={draw=white, fill=white}"))
	push!(G, Axis(plots, title="divided rectangles lower bound"))
end
plot(p)
# -

p = let
	function _get_shubert_piyavskii_intersection(A, B, L)
	    t = ((A[2] - B[2]) - L*(A[1] - B[1])) / (2L)
	    (A[1] + t, A[2] - t*L)
	end

	f = x -> log(1 + exp(x))/3.5 + exp(-(x+3)^2/16) + 3
	a, b = -5.0, 5.0

	# x_arr = [-3.2, 0, 1.5, 4.5]
	x_arr = [-3.33333, -1.11111, 0.0, 1.11111, 3.33333]
	y_arr = f.(x_arr)

	plots = Plots.Plot[]
	push!(plots, Plots.Linear(f, (a, b), style="black, thick, solid, mark=none", legendentry=L"f(x)"))

	L_arr = range(0.1,stop=1.0,length=10)
	for (iL, L) in enumerate(L_arr)
	    lowerbound_style = "solid, thick, pastelBlue!$(108-8iL), mark=none"
	    x_min, y_min = a, y_arr[1] - L*abs(a-x_arr[1])
	    push!(plots, Plots.Linear([a, x_arr[1]], [y_min, y_arr[1]], style=lowerbound_style, legendentry="\$\\ell = $L\$"))
	    for i in 2 : length(x_arr)
	        A = (x_arr[i-1], y_arr[i-1])
	        B = (x_arr[i],   y_arr[i])
	        if L*(B[1] - A[1]) > abs(B[2] - A[2])
	            C = _get_shubert_piyavskii_intersection(A, B, L)
	            if C[2] < y_min
	                x_min, y_min = C
	            end
	            push!(plots, Plots.Linear([A[1], C[1], B[1]], [A[2], C[2], B[2]], style=lowerbound_style * ", forget plot"))
	        end
	    end
	    push!(plots, Plots.Linear([x_arr[end], b], [y_arr[end], y_arr[end] - L*abs(b-x_arr[end])], style=lowerbound_style * ", forget plot"))
	    push!(plots, Plots.Scatter([x_min], [y_min], style="only marks, solid, mark=*, mark size=1.25, mark options={fill=pastelRed, draw=pastelRed}, forget plot"))
	    if iL == length(L_arr)
	        plots[end].style = plots[end].style[1:end-13]
	        plots[end].legendentry="lower-bound minima"
	    end
	end

	push!(plots, Plots.Scatter(x_arr, y_arr, style="only marks, solid, black, mark=*, mark size=1.25, mark options={fill=black, draw=black}", legendentry="evaluation points"))

	Axis(plots, xlabel=L"x", ylabel=L"y", xmin=a, xmax=b,
		 width="10cm", height="6cm",
		 style="legend pos=outer north east, clip marker paths, axis on top")
end
plot(p)

# +

struct Cell
    x::Float64
    y::Float64
    w::Float64
end

p = let
	f = x -> log(1 + exp(x))/3.5 + exp(-(x+3)^2/16) + 3
	a, b = -5.0, 5.0

	function Base.split(cell::Cell)
	    x,y,w = cell.x,cell.y,cell.w
	    x_left = x - 2w/6
	    x_right = x + 2w/6
	    return (Cell(x_left, f(x_left), w/3),
	            Cell(x, y, w/3),
	            Cell(x_right, f(x_right), w/3)
	           )
	end
	function split!(cells::Vector{Cell}, i::Int)
	    cell_a, cell_b, cell_c = split(cells[i])
	    insert!(cells, i, cell_a)
	    insert!(cells, i+2, cell_c)
	    cells[i+1] = cell_b
	    return cells
	end

	cells = [Cell(0.0, f(0.0), b-a)]
	split!(cells, 1)
	split!(cells, 2)

	plots = Plots.Plot[]
	push!(plots, Plots.Linear(f, (a, b), style="black, thick, solid, mark=none", legendentry=L"f(x)"))

	L_arr = range(0.1,stop=1.0,length=10)
	for (iL, L) in enumerate(L_arr)
	    lowerbound_style = "solid, thick, pastelBlue!$(108-8iL), mark=none"
	    x_min, y_min = NaN, Inf

	    for cell in cells
	        xl = max(a, cell.x - cell.w/2)
	        xr = min(b, cell.x + cell.w/2)
	        yl = cell.y - L*abs(cell.x - xl)
	        yr = cell.y - L*abs(cell.x - xr)

	        if yl < y_min
	            x_min, y_min = xl, yl
	        elseif yr < y_min
	            x_min, y_min = xr, yr
	        end

	        push!(plots, Plots.Linear([xl, cell.x, xr], [yl, cell.y, yr], style=lowerbound_style * ", forget plot"))
	    end
	    plots[end].style = plots[end].style[1:end-13]
	    plots[end].legendentry="\$\\ell = $L\$"

	    push!(plots, Plots.Scatter([x_min], [y_min], style="only marks, solid, mark=*, mark size=1.25, mark options={fill=pastelRed, draw=pastelRed}, forget plot"))
	    if iL == length(L_arr)
	        plots[end].style = plots[end].style[1:end-13]
	        plots[end].legendentry="lower-bound minima"
	    end
	end

	push!(plots, Plots.Scatter([cell.x for cell in cells], [cell.y for cell in cells], style="only marks, solid, black, mark=*, mark size=1.25, mark options={fill=black, draw=black}", legendentry="evaluation points"))

	Axis(plots, xlabel=L"x", ylabel=L"y", xmin=a, xmax=b,
		  width="10cm", height="6cm",
		 style="legend pos=outer north east, clip marker paths, axis on top")
end
plot(p)

# +
using PGFPlots

p = let

	plots = Plots.Plot[]
	x = 1/3
	y = 4.0
	L = 5.0
	push!(plots, Plots.Scatter([1,1,1/3,1/3,1/3,1/3,1/9,1/9,1/9,1/27,1/27,1/27],
	                           [10,9,8.5,5,4.7,4,6.2,4.1,3.8,4.15,3.7,3.2], style="mark=*, mark size=2, mark options={draw=black, fill=black}"))
	push!(plots, Plots.Linear([0,1.3], [y-L*x,y+L*(1.3-x)], style="black, solid, mark=none"))
	push!(plots, Plots.Linear([0,0.5], [y-L*x,y-L*x], style="gray, solid, mark=none"))
	push!(plots, Plots.Node(L"f(c_i) - \ell(b_i-a_i)/2",0.5, y-L*x, style="right"))
	push!(plots, Plots.Linear([0,0.5], [y,y], style="pastelBlue, solid, mark=none"))
	push!(plots, Plots.Node(L"f(c_i)",0.5, y, style="right"))
	push!(plots, Plots.Linear([x,x], [0,10.5], style="pastelBlue, solid, mark=none"))
	push!(plots, Plots.Node(L"(b_i-a_i)/2",0.4,0.6, axis="axis description cs"))
	push!(plots, Plots.Node("slope = \$\\ell\$",0.7,0.5, axis="axis description cs"))

	Axis(plots, width="8cm", ymin=1, xlabel="interval half-width \$(b-a)/2\$", ylabel="interval center value \$f(c)\$", style="xtick=\\empty, ytick=\\empty, axis lines=left")
end

plot(p)

# +
using PGFPlots

p = let

	plots = Plots.Plot[]
	x = 1/3
	y = 4.0
	L = 5.0
	push!(plots, Plots.Scatter([1,1,1/3,1/3,1/3,1/3,1/9,1/9,1/9,1/27,1/27,1/27],
	                           [10,9,8.5,5,4.7,4,6.2,4.1,3.8,4.15,3.7,3.2], style="mark=*, mark size=2, mark options={only marks, draw=black, fill=black}", legendentry="suboptimal"))
	push!(plots, Plots.Linear([1,1/3,1/27],
	                           [9,4,3.2], style="solid, pastelBlue, mark=*, mark size=2, mark options={draw=pastelBlue, fill=pastelBlue}", legendentry="potentially optimal"))

	Axis(plots, width="8cm", xmin=0, ymin=1, xlabel="interval half-width \$(b-a)/2\$", ylabel="interval center value \$f(c)\$", style="xtick=\\empty, ytick=\\empty, axis lines=left, legend cell align=left, legend style={at={(axis description cs:0.6,0.25)},anchor=west}")
end

plot(p)

# +
using DataStructures
struct Interval
    c
    y
    depths
end
const Intervals = OrderedDict{Float64,PriorityQueue{Interval, Float64}}

p = let

	f = x -> sin(x) + sin(2x) + sin(4x) + sin(8x)

	basis(i, n) = [k == i ? 1. : 0. for k in 1 : n]

	rev_unit_hypercube_parameterization(x, a, b) = x.*(b-a) + a
	function reparameterize_to_unit_hypercube(f, a, b)
	    Δ = b-a
	    return x->f(x.*Δ + a)
	end

	min_depth(interval) = minimum(interval.depths)
	vertex_dist(interval) = norm(0.5*3.0.^(-interval.depths), 2)

	function add_interval!(intervals, interval)
		d = min_depth(interval)
	    if !haskey(intervals, d)
	        intervals[d] = PriorityQueue{Interval, Float64}()
	    end
	    enqueue!(intervals[d], interval, interval.y)
	end

	function get_opt_intervals(intervals, ϵ, y_best)
	    stack = Interval[]
	    for (d, pq) in intervals
	    	if !isempty(pq)
	            interval = DataStructures.peek(pq)[1]
	            x, y = d, interval.y
	            while !isempty(stack)
	            	interval1 = stack[end]
	            	x1, y1 = vertex_dist(interval1), interval1.y
	            	L1 = (y - y1)/(x - x1)
	            	if y1 - L1*x1 > y_best - ϵ || y < y1
	            		pop!(stack)
	            	elseif length(stack) > 1
	            		interval2 = stack[end-1]
	            		x2, y2 = vertex_dist(interval2), interval2.y
	            		L2 = (y1 - y2)/(x1 - x2)
	            		if L2 > L1
	            			pop!(stack)
	                    else
	                        break
	            		end
	                else
	                    break
	            	end
	           	end
	            push!(stack, interval) # add new interval
	        end
	    end
	    stack
	end

	function divide(f, interval)
	    c, d, n = interval.c, min_depth(interval), length(interval.c)
	    dirs = findall(interval.depths .== d)
	    cs = [(c + 3.0^(-d-1)*basis(i,n),
	           c - 3.0^(-d-1)*basis(i,n)) for i in dirs]
	    vs = [(f(C[1]), f(C[2])) for C in cs]
	    minvals = [min(V[1], V[2]) for V in vs]

	    retval = Interval[]
	    depths = copy(interval.depths)
	    for j in sortperm(minvals)
	        depths[dirs[j]] += 1
	        C, V = cs[j], vs[j]
	        push!(retval, Interval(C[1], V[1], copy(depths)))
	        push!(retval, Interval(C[2], V[2], copy(depths)))
	    end
	    push!(retval, Interval(c, interval.y, copy(depths)))
	    retval
	end

	function add_axes!(g, intervals, opt_intervals)

	    x_pts = collect(range(-2,stop=2,length=200))
	    y_pts = f.(x_pts)

	    style_normal = "clip marker paths, solid, only marks, mark=*, mark size=1, mark options={draw=black, fill=black}"
	    style_opt = "clip marker paths, solid, only marks, mark=*, mark size=1.25, mark options={draw=pastelBlue, fill=pastelBlue}"

	    arr_x_normal = Float64[]
	    arr_y_normal = Float64[]
	    arr_h_normal = Float64[]
	    arr_x_opt = Float64[]
	    arr_y_opt = Float64[]
	    arr_h_opt = Float64[]

	    plot_intervals = Plots.Plot[]
	    push!(plot_intervals, Plots.Linear(x_pts, y_pts, style="black, solid, mark=none"))
	    for (d,Is) in intervals
	        for I in keys(Is)
	            x = rev_unit_hypercube_parameterization(I.c, a, b)[1]
	            y = I.y
	            h = 2*3.0^(-d)
	            if in(I, opt_intervals)
	                push!(arr_x_opt, x)
	                push!(arr_y_opt, y)
	                push!(arr_h_opt, h)
	                push!(plot_intervals, Plots.Linear([x-h,x+h],[y,y], style="solid, pastelBlue, mark=none"))
	            else
	                push!(arr_x_normal, x)
	                push!(arr_y_normal, y)
	                push!(arr_h_normal, h)
	                push!(plot_intervals, Plots.Linear([x-h,x+h],[y,y], style="solid, gray, mark=none"))
	            end
	        end
	    end

	    push!(plot_intervals, Plots.Scatter(arr_x_normal, arr_y_normal, style=style_normal))
	    push!(plot_intervals, Plots.Scatter(arr_x_opt, arr_y_opt, style=style_opt))

	    push!(g, Axis(plot_intervals, xmin=-2, xmax=2, ymin=-2.5, ymax=2.5, ylabel=L"y"))

	    p_scatter_normal = Plots.Scatter(arr_h_normal, arr_y_normal, style=style_normal)
	    p_scatter_opt = Plots.Scatter(arr_h_opt, arr_y_opt, style=style_opt)

	    push!(g, Axis([p_scatter_normal, p_scatter_opt], xmin=0.0, xmax=2.0, ymin=-2.5, ymax=2.5))

	    g
	end

	ϵ = 0.001
	a = [-2.0]
	b = [2.0]
	g = reparameterize_to_unit_hypercube(x->f(x[1]), a, b)
	intervals = Intervals()
	n = length(a)
	c = fill(0.5, n)
	interval = Interval(c, g(c), fill(0, n))
	add_interval!(intervals, interval)
	c_best, y_best = copy(interval.c), interval.y

	K = 6
	G = GroupPlot(2, K, groupStyle="horizontal sep=0.25cm, vertical sep=0.25cm, xlabels at=edge bottom, xticklabels at=edge bottom, ylabels at=edge left, yticklabels at=edge left", style="width=6cm, height=3.25cm")

	for k in 1 : K
	    S = get_opt_intervals(intervals, ϵ, y_best)

	    add_axes!(G, intervals, get_opt_intervals(intervals, ϵ, y_best))

	    to_add = Interval[]
	    for interval in S
	        append!(to_add, divide(g, interval))
	        dequeue!(intervals[min_depth(interval)])
	    end
	    for interval in to_add
	        add_interval!(intervals, interval)
	        if interval.y < y_best
	            c_best, y_best = copy(interval.c), interval.y
	        end
	    end
	end

	G.axes[end-1].xlabel=L"x"
	G.axes[end].xlabel="interval half-width"
	G
end

plot(p)
# -

	using PGFPlots
	using DataStructures
	using Colors

	p = let

		function reparameterize_to_unit_hypercube(f, a, b)
		    Δ = b-a
		    x->f(x.*Δ + a)
		end

		min_depth(interval) = minimum(interval.depths)
		vertex_dist(interval) = norm(0.5*3.0.^(-interval.depths), 2)

		function add_interval!(intervals, interval)
			d = vertex_dist(interval)
		    if !haskey(intervals, d)
		        intervals[d] = PriorityQueue{Interval, Float64}()
		    end
		    enqueue!(intervals[d], interval, interval.y)
		end

		function get_opt_intervals(intervals, ϵ, y_best)
		    stack = Interval[]
		    for (x, pq) in intervals
		    	if !isempty(pq)
		            interval = DataStructures.peek(pq)[1]
		            y = interval.y

		            # If no previous point, we keep it.
		            # If one previous point, keep it if new one is lower.
		            # If more than one previous point
		            #    - reject if new point is any higher than prev point
		            #    - keep otherwise
		            #        while has two prev points:
		            #        - reject prev point if it falls below line with new and its prev pt.

		            while length(stack) > 1
		                interval1 = stack[end]
		                interval2 = stack[end-1]
		                x1, y1 = vertex_dist(interval1), interval1.y
		                x2, y2 = vertex_dist(interval2), interval2.y
		                ℓ = (y2 - y) / (x2 - x)
		                if y1 <= ℓ*(x1-x) + y + ϵ
		                    break
		                end
		                # Remove previous interval
		                pop!(stack)
		            end

		            if !isempty(stack) && interval.y > stack[end].y + ϵ
		                # Skip new interval
		                continue
		            end

		            push!(stack, interval) # add new interval
		        end
		    end
		    return stack
		end

		function basis(i, n)
		    eᵢ = zeros(n)
		    eᵢ[i] = 1
		    eᵢ
		end
		function divide(f, interval)
		    c, d = interval.c, min_depth(interval)
		    n = length(c)
		    dirs = findall(interval.depths .== d)
		    cs = [(c + 3.0^(-d-1)*basis(i,n),
		           c - 3.0^(-d-1)*basis(i,n)) for i in dirs]
		    vs = [(f(C[1]), f(C[2])) for C in cs]
		    minvals = [min(V[1], V[2]) for V in vs]

		    retval = Interval[]
		    depths = copy(interval.depths)
		    for j in sortperm(minvals)
		        depths[dirs[j]] += 1
		        C, V = cs[j], vs[j]
		        push!(retval, Interval(C[1], V[1], copy(depths)))
		        push!(retval, Interval(C[2], V[2], copy(depths)))
		    end
		    push!(retval, Interval(c, interval.y, copy(depths)))

		    retval
		end

		function plot_topdown(rects, g)
		    arr_x = Float64[]
		    arr_y = Float64[]
		    for list in values(rects)
		        for r in keys(list)
		            push!(arr_x, r.c[1])
		            push!(arr_y, r.c[2])
		        end
		    end

		    p = Plots.Plot[]
		    push!(p, PGFPlots.Plots.Image((x,y)->g([x,y]), (0,1), (0,1), colormap=pasteljet, xbins=600, ybins=600))
        #   push!(p, PGFPlots.Plots.Contour(g, (0,1), (0,1), levels=[0.1,0.2,0.5,1,2,3,4]))
		#     push!(p, Plots.Scatter(arr_x, arr_y, style="mark=*, mark size=1, mark options={draw=black, fill=black}"))

		    for list in values(rects)
		        for r in keys(list)
		            x_lo = r.c[1] - 0.5*3.0^(-r.depths[1])
		            x_hi = r.c[1] + 0.5*3.0^(-r.depths[1])
		            y_lo = r.c[2] - 0.5*3.0^(-r.depths[2])
		            y_hi = r.c[2] + 0.5*3.0^(-r.depths[2])
		            push!(p, Plots.Linear([x_lo,x_lo,x_hi,x_hi,x_lo],[y_lo,y_hi,y_hi,y_lo,y_lo], style="solid, white, mark=none"))
		        end
		    end

		    Axis(p, xlabel=L"x_1", ylabel=L"x_2", axisEqual=true, width="10cm", height="10cm", enlargelimits=false,
		            style="xtick={0,1}, ytick={0,1}")
		end

		f = x-> (-1.275*x[1]^2/(π^2) + 5*x[1]/π + x[2] - 6)^2 + (10 - 5/(4π))*cos(x[1]) + 10
		g = reparameterize_to_unit_hypercube(f, [-5.0,0.0], [10.0,15.0])

		intervals = OrderedDict{Float64, PriorityQueue{Interval, Float64}}()
		ϵ = 0.01

		c = [0.5,0.5]
		n = length(c)
		interval = Interval(c, g(c), fill(0, n))
		add_interval!(intervals, interval)
		y_best = interval.y

		for k in 1 : 16
		    S = get_opt_intervals(intervals, ϵ, y_best)
		    to_add = Interval[]
		    for r in S
		        append!(to_add, divide(g, r))
		        dequeue!(intervals[vertex_dist(r)])
		    end
		    for r in to_add
		        add_interval!(intervals, r)
		        if interval.y < y_best
	                y_best = interval.y
	            end
		    end
		end

		plot_topdown(intervals, g)
	end

	plot(p)

p = let
	plots = Plots.Plot[]
	push!(plots, Plots.Scatter([1/2,1/6,1/6],
	                           [1.231,2.029,1.861], style="only marks, mark=*, mark size=2, mark options={draw=black, fill=black}", legendentry="suboptimal"))
	push!(plots, Plots.Linear([1/2,1/6],
                               [0.5,0.158], style="solid, pastelBlue, mark=*, mark size=2, mark options={draw=pastelBlue, fill=pastelBlue}", legendentry="potentially optimal"))
    push!(plots, Plots.Node(L"1",0.5,0.5, style="right"))
    push!(plots, Plots.Node(L"2",0.5,1.231, style="right"))
    push!(plots, Plots.Node(L"3",1/6,0.158, style="left"))
    push!(plots, Plots.Node(L"4",1/6,2.029, style="left, yshift= 2pt"))
    push!(plots, Plots.Node(L"5",1/6,1.861, style="left, yshift=-2pt"))


	Axis(plots, width="6cm", xmin=0, xmax=0.6, ymin=0, ymax=2.5, xlabel="interval half-width", ylabel="interval center value", style="axis lines=left, legend cell align=left, xtick={0,0.166666,0.5}, xticklabels={\$0\$,\$1/6\$,\$1/2\$}, legend style={at={(axis description cs:1.25,0.25)},anchor=west}")
end
plot(p)

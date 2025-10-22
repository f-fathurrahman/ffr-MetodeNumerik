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

# # Discrete
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

        using LinearAlgebra
        p = let
            f_ = (x,y) -> x + y
            dom = (-3.25,3.25)

            Xgood = Vector{Float64}[]
            Xbad = Vector{Float64}[]
            for x in -3 : 3
                for y in -3 : 3
                    p = [x,y]
                    if norm(p) ≤ 2
                        push!(Xgood, p)
                    else
                        push!(Xbad, p)
                    end
                end
            end

            y_arr = (x->f_(x[1],x[2])).(Xgood)
            y_min = minimum(y_arr)
            Xopt = Xgood[findall(y_arr .== y_min)]

            Axis([
                Plots.Image(f_, dom, dom, xbins=600, ybins=600, colormap=pasteljet, colorbar = false),
                Plots.Contour((x,y)->norm([x,y]), dom, dom, levels=[2], contour_style="draw color=white", style="thick, white", xbins=101, ybins=101),
                Plots.Scatter([x[1] for x in Xgood], [x[2] for x in Xgood], style="solid, white, mark=*, mark size=1.25, mark options={draw=white, fill=white}"),
                Plots.Scatter([x[1] for x in Xbad],  [x[2] for x in Xbad], style="solid, white, mark=o, mark size=1.25, mark options={draw=white}"),
                Plots.Scatter([x[1] for x in Xopt],  [x[2] for x in Xopt], style="solid, pastelRed, mark=*, mark size=1.5, mark options={draw=pastelRed, fill=pastelRed}"),
                Plots.Scatter([-sqrt(2)],  [-sqrt(2)], style="solid, pastelRed, mark=+, mark size=5, mark options={draw=pastelRed, fill=pastelRed}"),
            ], xlabel=L"x_1", ylabel=L"x_2", width="4.75cm", height="4.75cm", style="contour/labels=false, view={0}{90}")
        end
        plot(p)

# +
p = let
f_ = x->1 - exp(-2(x-2)^2) - 1.25exp(-2(x-4)^2)
dom = (0,5)
m = 11
x_pts = collect(range(dom[1], stop=dom[2], length=m))

a = 0.0
b = 3.75
lo = -0.5
hi =  1.25

Axis([
Plots.Linear(f_, dom, style="solid, black, mark=none"),
Plots.Linear([a,a], [lo,hi], style="name path=A, draw=none, mark=none"),
Plots.Linear([b,b], [lo,hi], style="name path=B, draw=none, mark=none"),
Plots.Command("\\addplot[pastelBlue!40] fill between[of=A and B]"),
Plots.Scatter(x_pts[1:8], f_.(x_pts[1:8]), style="solid, black, mark=*, mark size=1, mark options={draw=black, fill=black}"),
Plots.Scatter(x_pts[9:end], f_.(x_pts[9:end]), style="solid, black, mark=*, mark size=1.25, mark options={draw=black, fill=white}"),
Plots.Scatter([b], [f_(b)], style="solid, pastelRed, mark=*, mark size=1.25, mark options={draw=pastelRed, fill=pastelRed}"),
# Plots.Scatter([2], [f_(2)], style="solid, white, mark=*, mark size=1.25, mark options={draw=white, fill=white}"),
Plots.Node(L"x^*", 2, f_(2), style="below"),
Plots.Command("\\draw[->, pastelRed] (axis cs:$b,$(f_(b))) -- (axis cs:3.5,$(f_(3.5)));"),
], style="enlarge x limits=0, xtick=\\empty, ytick=\\empty, xlabel=\$x\$, ylabel=\$y\$", ymin=lo, ymax=hi, width="9cm"
)
end
plot(p)
# -

    p = let
    f_ = x->sin(x)*sin(2x+0.1)*sin(4x+0.2)*sin(8x+0.3)
    dom = (0,3)
    m = 20
    x_pts = collect(range(dom[1],stop=dom[2],length=m))

    G = GroupPlot(2,1, groupStyle="horizontal sep=0.5cm, ylabels at=edge left",
    style="width=7cm, xtick=\\empty, ytick=\\empty, xlabel=\$x\$, ylabel=\$y\$")
    push!(G, Axis(
    Plots.Linear(f_, dom, style="solid, black, mark=none", xbins=201),
    style="enlarge x limits=0"
    ))
    push!(G, Axis(
    Plots.Scatter(x_pts, f_.(x_pts), style="solid, black, mark=*, mark size=1, mark options={draw=black, fill=black}")
    ))
    G
    end
    plot(p)

        using LightGraphs
        using StatsBase
        using Vec
        using Random

        function execute_ant(graph, lengths, A, τ, best_tour, best_tour_length)
            path = [1]
            while length(path) < nv(graph)
                i = path[end]
                neighbors = outneighbors(graph, i)
                valid_neighbors = setdiff(neighbors, path)
                if isempty(valid_neighbors) # ant got stuck
                    return (best_tour, best_tour_length)
                end
                as = [A[(i,j)] for j in valid_neighbors]
                push!(path, valid_neighbors[sample(Weights(as))])
            end

            l = sum(lengths[(path[i-1],path[i])] for i in 2 : length(path))
            for i in 2 : length(path)
                τ[(path[i-1],path[i])] += 1/l
            end
            if l < best_tour_length
                return (path, l)
            else
                return (best_tour, best_tour_length)
            end
        end

        function edge_attractiveness(graph, τ, η, α, β)
            A = Dict{Tuple{Int,Int},Float64}()
            for i in 1 : nv(graph)
                neighbors = outneighbors(graph, i)
                for j in neighbors
                    v = τ[(i,j)]^α * η[(i,j)]^β
                    A[(i,j)] = v
                end
            end
            return A
        end

        let
            dist_threshold = 0.8

            function closest_point(P::VecE2, vertices::Vector{VecE2})
                d_min = Inf
                retval = 0
                for (i,Q) in enumerate(vertices)
                    d = abs(Q-P)
                    if 0 < d < d_min
                        d_min = d
                        retval = i
                    end
                end
                return retval
            end
            function edge_collides_with_other_edges(i::Int, j::Int, vertices::Vector{VecE2}, graph::DiGraph)
                seg = LineSegment(vertices[i], vertices[j])
                for e in edges(graph)
                    A = vertices[e.src]
                    B = vertices[e.dst]
                    r = polar(1.0, atan(B - A))
                    seg2 = LineSegment(A + r*0.05, B - r*0.05)
                    if intersects(seg, seg2)
                        return true
                    end
                end
                return false
            end

            Random.seed!(1)
            vertices = VecE2[]
            for θ in range(0, stop=2π, length=21)[1:20]
                push!(vertices, polar(1.25 + sin(4θ),θ))
            end

            graph = DiGraph(length(vertices))
            for i in 1 : length(vertices)-1
                add_edge!(graph, i, i+1)
            end
            add_edge!(graph, length(vertices), 1)

            for k in 1 : 100
                i = rand(1:nv(graph))
                j = sample([j for j in 1 : nv(graph) if norm(vertices[i] - vertices[j]) < 2])
                if i != j && !has_edge(graph, i, j) && !has_edge(graph, j, i) &&
                             !edge_collides_with_other_edges(i, j, vertices, graph)

                    add_edge!(graph,i,j)
                end
            end
            for k in 1 : 14
                i = rand(1:nv(graph))
                j = sample([j for j in 1 : nv(graph)])
                if i != j && !has_edge(graph, i, j) && !has_edge(graph, j, i) &&
                             !edge_collides_with_other_edges(i, j, vertices, graph)

                    add_edge!(graph,i,j)
                end
            end

            lengths = Dict{Tuple{Int,Int}, Float64}()
            for e in edges(graph)
                lengths[(e.src, e.dst)] = norm(vertices[e.src] - vertices[e.dst])
            end

            Random.seed!(0)

            k_max = 1*4
            m = 50
            α = 1.0
            β = 5.0
            ρ = 0.5
            η = Dict((e.src,e.dst)=>1/lengths[(e.src,e.dst)] for e in edges(graph))
            τ = Dict((e.src,e.dst)=>1.0 for e in edges(graph))

            # rendering params
            r = 0.03
            scale = 0.7

            best_tour = Int[]
            best_tour_length = Inf

            Δx = 4.5
            Δy = 4.5
            max_node_x = 3
            nodex = 0
            nodey = 0

            o = stdout # IOBuffer()
            println(o, "\\begin{tikzpicture}[x=$(scale)cm, y=$(scale)cm]")

            for k in 1 : k_max

                println(o, "\\node at ($(Δx*nodex),$(Δy*nodey)) { \\begin{tikzpicture}")
                τmax = maximum(values(τ))
                for e in edges(graph)
                    A = vertices[e.src]
                    B = vertices[e.dst]
                    α = τ[(e.src, e.dst)] / τmax
                    println(o, "\\draw[-{Latex[length=0.7mm, width=0.5mm]}, shorten >=3pt, pastelBlue!$(round(90*α + 10))] ($(A.x), $(A.y)) -- ($(B.x), $(B.y));")
                end
                for A in vertices
                    println(o, "\\fill ($(A.x), $(A.y)) circle ($r);")
                end
                println(o, "\\end{tikzpicture} };")
                nodex += 1
                if nodex > max_node_x
                    nodex = 0
                    nodey -= 1
                end

                A = edge_attractiveness(graph, τ, η, α, β)
                for (e,v) in τ
                    τ[e] = (1-ρ)*v
                end
                for ant in 1 : m
                    (best_tour, best_tour_length) = execute_ant(graph, lengths, A, τ, best_tour, best_tour_length)
                end
            end

            println(o, "\\end{tikzpicture}")
        end

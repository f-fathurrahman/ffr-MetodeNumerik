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

# # Linear
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
using Vec
using LinearAlgebra

p = let

    g = GroupPlot(4,1,groupStyle="horizontal sep=0.5cm, ylabels at=edge left",
    	style = "width=5cm, height=5cm, xmin=0, xmax=1, ymin=0, ymax=1, ytick=\\empty, xtick=\\empty, xlabel = \$x_1\$, ylabel= \$x_2\$, contour/labels=false, axis equal, view={0}{90}"
    )



    c = Float64[2, 1]
    f = x->c⋅x

    A = VecE2(0.20,0.20)
    B = A + polar(2.0, deg2rad(10))
    C = A + polar(2.0, deg2rad(75))

    arrow_lo = [0.7,0.9]
    arrow_hi = arrow_lo - 0.6*normalize(c)
    push!(g, Axis([Plots.Contour(f, (0,1), (0,1), contour_style="draw color=gray"),
    			   Plots.Linear([A.x, B.x], [A.y, B.y], style="solid, black, mark=none, name path=A"),
    			   Plots.Linear([A.x, C.x], [A.y, C.y], style="solid, black, mark=none, name path=B"),
    			   Plots.Command("\\addplot[pastelBlue!50] fill between[of=A and B]"),
    			   Plots.Scatter([A.x], [A.y], style="mark=*, mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}"),
    			   Plots.Node(L"\vect x^*",A.x, A.y, style="below right"),
    			   Plots.Command("\\draw[myarrow] (axis cs:$(arrow_lo[1]), $(arrow_lo[2])) -- (axis cs:$(arrow_hi[1]), $(arrow_hi[2]));"),
    			   Plots.Node(L"-\vect c", (arrow_lo[1] + arrow_hi[1])/2, (arrow_lo[2] + arrow_hi[2])/2, style="below right", axis="axis description cs"),
    		], title="One Solution"
    ))

    r = VecE2(1,-2)
    q = VecE2(2, 1)

    A = VecE2(0.20,0.20)
    B = A + polar(2.0, deg2rad(10))
    C = A + polar(2.0, deg2rad(75))
    push!(g, Axis([Plots.Contour(x->-f(x), (0,1), (0,1), contour_style="draw color=gray"),
    			   Plots.Linear([A.x, B.x], [A.y, B.y], style="solid, black, mark=none, name path=A"),
    			   Plots.Linear([A.x, C.x], [A.y, C.y], style="solid, black, mark=none, name path=B"),
    			   Plots.Command("\\addplot[pastelBlue!50] fill between[of=A and B]"),
    			   Plots.Command("\\draw[myarrow] (axis cs:$(arrow_hi[1]), $(arrow_hi[2])) -- (axis cs:$(arrow_lo[1]), $(arrow_lo[2]));"),
    			   Plots.Node(L"-\vect c", (arrow_lo[1] + arrow_hi[1])/2, (arrow_lo[2] + arrow_hi[2])/2, style="below right", axis="axis description cs"),
    		], title="Unbounded Solution"
    ))

    A = 0.1*q - 0.15*r
    B = 0.1*q + 0.00*r
    C = B + polar(0.75, deg2rad(10))
    D = A + polar(0.6, deg2rad(75))
    push!(g, Axis([Plots.Contour(f, (0,1), (0,1), contour_style="draw color=gray"),
                   Plots.Linear([A.x, B.x, C.x, D.x, A.x],[A.y, B.y, C.y, D.y, A.y], style="solid, black, mark=none, name path=A"),
                   Plots.Linear([A.x, B.x],[A.y, B.y], style="solid, ultra thick, pastelBlue, mark=none, name path=B"),
                   Plots.Command("\\addplot[pastelBlue!50] fill between[of=A and B]"),
    			   Plots.Node(L"\vect x^*",(A.x+B.x)/2, (A.y+B.y)/2, style="above right"),
    			   Plots.Command("\\draw[myarrow] (axis cs:$(arrow_lo[1]), $(arrow_lo[2])) -- (axis cs:$(arrow_hi[1]), $(arrow_hi[2]));"),
    			   Plots.Node(L"-\vect c", (arrow_lo[1] + arrow_hi[1])/2, (arrow_lo[2] + arrow_hi[2])/2 + 0.05, style="below right", axis="axis description cs"),
    		], title="Infinite Solutions"
    ))

    push!(g, Axis([Plots.Contour(f, (0,1), (0,1), contour_style="draw color=gray"),
    			   Plots.Command("\\draw[myarrow] (axis cs:$(arrow_lo[1]), $(arrow_lo[2])) -- (axis cs:$(arrow_hi[1]), $(arrow_hi[2]));"),
    			   Plots.Node(L"-\vect c", (arrow_lo[1] + arrow_hi[1])/2, (arrow_lo[2] + arrow_hi[2])/2, style="below right", axis="axis description cs"),
    		], title="No Solution"
    ))

    g
end

plot(p)
# -

p = let
	Axis(
		[Plots.Linear([0,1], [-1,0], style="solid, black, mark=none"),
		 Plots.Linear([1,4], [0,3], style="solid, thick, pastelBlue, mark=none"),
		 ],
		xlabel=L"x", ylabel=L"s",
		xmin=0, xmax=4.5, ymin=-1.25, ymax=3.25, width="1.2*9cm",
		style="axis lines=center",
	)
end
plot(p)

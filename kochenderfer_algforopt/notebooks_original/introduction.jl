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

# # Introduction
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +

p = let
	f = (x,y) -> (1-x)^2 + 5*(y - x^2)^2
	Axis([Plots.Contour(f, (-2,2), (-2,2), levels=[1,2,3,5,10,20,50,100]),
		  Plots.Scatter([1], [1], style="black, solid, mark=*, mark size = 1, mark options={draw=black}"),
		  ], width="9cm", height="9cm", xlabel=L"x_1", ylabel=L"x_2", style="view={0}{90}, contour/labels=false"
	)
end

plot(p)
# -

	p = let
		xdom = (-2.5,2.5)
		ydom = (-2.5,2.5)

		g = GroupPlot(2,1,groupStyle="horizontal sep=1.5cm, ylabels at=edge left",
                  	  style="width=5.25cm, xlabel=\$x_1\$, ylabel=\$x_2\$, xmin=$(xdom[1]), xmax=$(xdom[2]), ymin=$(ydom[1]), ymax=$(ydom[2])")

		f = (x,y) -> x^2 - y^2

		c1 = x -> -sqrt(2 + x^2)
		c2 = x -> -sqrt(4 + x^2)
		c3 = y ->  sqrt(2 + y^2)
		c4 = y ->  sqrt(4 + y^2)

		xs = range(xdom[1],stop=xdom[2],length=51)
		xs3 = range(-2, stop=0.5, length=41)
		xs4 = range(-1.5, stop=0, length=41)

		push!(g,
			Axis(
				[
					Plots.Command("""
					    \\addplot3[patch,patch refines=3,
							shader=faceted interp,
							patch type=biquadratic]
					    	table[z expr=x^2-y^2]
					    {
					        x  y
					        -2.5 -2.5
					        2.5  -2.5
					        2.5  2.5
					        -2.5 2.5
					        0  -2.5
					        2.5  0
					        0  2.5
					        -2.5 0
					        0  0
					    };
					"""),
					Plots.Linear3([-2.5,0.7],[-2.5,0.7],[0,0], style="solid, white, mark=none"),
					Plots.Linear3([-1.5,2.5],[1.5,-2.5],[0,0], style="solid, white, mark=none"),
					Plots.Linear3(xs,[c1(x) for x in xs],fill(-2, length(xs)), style="solid, white, mark=none"),
					Plots.Linear3(xs,[c2(x) for x in xs],fill(-4, length(xs)), style="solid, white, mark=none"),
					Plots.Linear3([ c3(x) for x in xs3], xs3, fill( 2, length(xs3)), style="solid, white, mark=none"),
					Plots.Linear3([ c4(x) for x in xs4], xs4, fill( 4, length(xs4)), style="solid, white, mark=none"),
					Plots.Linear3([-c3(x) for x in xs], xs, fill( 2, length(xs)), style="solid, white, mark=none"),
					Plots.Linear3([-c4(x) for x in xs], xs, fill( 4, length(xs)), style="solid, white, mark=none"),
				]
			)
		)
		# push!(g,
		# 	Axis(Plots.Image(f, xdom, ydom, xbins=600, ybins=600, colormap = pasteljet, colorbar = false),
		# 	     xmin=xdom[1], xmax=xdom[2], ymin=ydom[1], ymax=ydom[2], width="5.25cm", height="5.25cm", style="view={0}{90}", xlabel=L"x_1", ylabel=L"x_2")
		# )
		push!(g,
			Axis(Plots.Contour(f, xdom, ydom, xbins=100, ybins=100, levels=[-4,-2,0,2,4]),
			     style="axis equal, height=5.25cm, view={0}{90}")
		)

		g
	end

	plot(p)

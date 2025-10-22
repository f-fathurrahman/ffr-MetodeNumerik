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

# # Mdo
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

p = let

	g = GroupPlot(1, 2, groupStyle="vertical sep=1em, xlabels at=edge bottom, ylabels at=edge left, xticklabels at=edge bottom, yticklabels at=edge bottom", style="width=10cm, height=3.75cm, xlabel=iteration")

	function F1(A)
		A["y1"] = A["y2"] - A["x"]
		return A
	end
	function F2(A)
		A["y2"] = sin(A["y1"] + A["y3"])
		return A
	end
	function F3(A)
	    A["y3"] = cos(A["x"] + A["y2"] + A["y1"])
	    return A
	end

	_K = 20
	arr_t = collect(0:_K)

	function get_plot(Fs, add_legend)

		A = Dict{String,Float64}("x"=>1.0, "y1"=>1.0, "y2"=>1.0, "y3"=>1.0)

		arr_y1 = Array{Float64}(undef, _K+1)
		arr_y2 = Array{Float64}(undef, _K+1)
		arr_y3 = Array{Float64}(undef, _K+1)
		arr_y1[1] = A["y1"]
		arr_y2[1] = A["y2"]
		arr_y3[1] = A["y3"]

		for i in 1 : _K
			for F in Fs
				F(A)
			end
			arr_y1[i+1] = A["y1"]
			arr_y2[i+1] = A["y2"]
			arr_y3[i+1] = A["y3"]
		end
		p1 = Plots.Linear(arr_t, arr_y1, style="pastelPurple, solid, mark=*, mark size=1, mark options={draw=pastelPurple, fill=pastelPurple}")
		p2 = Plots.Linear(arr_t, arr_y2, style="pastelBlue, solid, mark=*, mark size=1, mark options={draw=pastelBlue, fill=pastelBlue}")
		p3 = Plots.Linear(arr_t, arr_y3, style="pastelSeaGreen, solid, mark=*, mark size=1, mark options={draw=pastelSeaGreen, fill=pastelSeaGreen}")
		if add_legend
			p1.legendentry = L"y^{(1)}"
			p2.legendentry = L"y^{(2)}"
			p3.legendentry = L"y^{(3)}"
		end


		return Axis([p1, p2, p3], style= add_legend ? "legend pos=outer north east" : nothing)
	end

	push!(g, get_plot([F1, F2, F3], true))
	push!(g, get_plot([F1, F3, F2], false))
	g
end
plot(p)

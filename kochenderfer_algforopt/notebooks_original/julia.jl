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

# # Julia
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

    		println("\\begin{tabular}{ll}")
    		println("\\toprule")
    		println("Package & Version \\\\")
    		println("\\midrule")
    		open("REQUIRE") do f
    		for l in readlines(f)
    		if startswith(l, "julia") || startswith(l, "Vec") || startswith(l, "Weave") || startswith(l, "ColorSchemes") || startswith(l, "Optim") || startswith(l, "PGFPlots") || startswith(l, "SymEngine") || startswith(l, "TikzGraphs") || length(l) < 1
    		continue
    		end
    		println("$l & $(Pkg.installed(l))\\\\")
    		end
    		end
    		println("\\bottomrule")
    		println("\\end{tabular}")

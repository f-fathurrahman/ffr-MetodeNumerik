
import QuanticsGrids as QG
import TensorCrossInterpolation as TCI
using QuanticsTCI: QuanticsTCI, quanticscrossinterpolate, integral

using Plots
using PlotThemes
themes!(:dark)

B = 2^(-30);
function f(x)
    return cos(x / B) * cos(x / (4 * sqrt(5) * B)) * exp(-x^2) + 2 * exp(-x)
end

xs = LinRange(0, 2.0^(-23), 1000)
plt = plot(title = "$(nameof(f))")
plot!(plt, xs, f.(xs), label = "$(nameof(f))", legend = true)


xs2 = LinRange(2.0^(-23), 3, 100000)
plt = plot(title = "$(nameof(f))")
plot!(plt, xs2, f.(xs2), label = "$(nameof(f))", legend = true)

R = 40 # number of bits
xmin = 0.0
xmax = log(20.0)
N = 2^R # size of the grid
# * Uniform grid (includeendpoint=false, default):
#   -xmin, -xmin+dx, ...., -xmin + (2^R-1)*dx
#     where dx = (xmax - xmin)/2^R.
#   Note that the grid does not include the end point xmin.
#
# * Uniform grid (includeendpoint=true):
#   -xmin, -xmin+dx, ...., xmin-dx, xmin,
#     where dx = (xmax - xmin)/(2^R-1).
qgrid = QG.DiscretizedGrid{1}(R, xmin, xmax; includeendpoint = true)
ci, ranks, errors = quanticscrossinterpolate(Float64, f, qgrid; maxbonddim = 15)


for i in [1, 2, 3, 2^R] # Linear indices
    # restore original coordinate `x` from linear index `i`
    x = QG.grididx_to_origcoord(qgrid, i)
    println("x: $(x), i: $(i), tci: $(ci(i)), ref: $(f(x))")
end

maxindex = QG.origcoord_to_grididx(qgrid, 2.0^(-23))
testindices = Int.(round.(LinRange(1, maxindex, 1000)))

xs = [QG.grididx_to_origcoord(qgrid, i) for i in testindices]
ys = f.(xs)
yci = ci.(testindices)

plt = plot(title = "$(nameof(f)) and TCI", xlabel = "x", ylabel = "y")
plot!(plt, xs, ys, label = "$(nameof(f))", legend = true)
plot!(plt, xs, yci, label = "tci", linestyle = :dash, alpha = 0.7, legend = true)

maxindex = QG.origcoord_to_grididx(qgrid, 2.0^(-23))
testindices = Int.(round.(LinRange(1, maxindex, 1000)))

xs = [QG.grididx_to_origcoord(qgrid, i) for i in testindices]
ys = f.(xs)
yci = ci.(testindices)
plt = plot(title = "x vs interpolation error: $(nameof(f))",
    xlabel = "x", ylabel = "interpolation error")
semilogy!(xs, abs.(ys .- yci), label = "log(|f(x) - ci(x)|)", yscale = :log10,
    legend = :bottomright, ylim = (1e-16, 1e-7), yticks = 10.0 .^ collect(-16:1:-7))


plt = plot(title = "x vs interpolation error: $(nameof(f))",
    xlabel = "x", ylabel = "interpolation error")

testindices = Int.(round.(LinRange(1, 2^R, 1000)))
xs = [QG.grididx_to_origcoord(qgrid, i) for i in testindices]
ys = f.(xs)
yci = ci.(testindices)
semilogy!(xs, abs.(ys .- yci), label = "log(|f(x) - ci(x)|)", legend = true,
    ylim = (1e-16, 1e-6), yticks = 10.0 .^ collect(-16:1:-6))


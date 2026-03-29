import MyQuanticsGrids as QG
import MyTensorCrossInterpolation as TCI
using MyQuanticsTCI: quanticscrossinterpolate, integral

B = 2^(-30);
function f(x)
    return cos(x / B) * cos(x / (4 * sqrt(5) * B)) * exp(-x^2) + 2 * exp(-x)
end

R = 40 # number of bits
xmin = 0.0
xmax = log(20.0)
N = 2^R # size of the grid
qgrid = QG.DiscretizedGrid{1}(R, xmin, xmax; includeendpoint = true)
ci, ranks, errors = quanticscrossinterpolate(Float64, f, qgrid; maxbonddim = 15)

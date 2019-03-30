# Numerical integration using trapezoidal rule
#
# f : function to be integrated, defined as a \texttt{Function}
# a : lower bound
# b : upper bound
# B : number
function integ_trapezoid( f::Function, a, b, N )
    x = a
    h = (b-a)/N
    s = f(x)
    #
    for i = 1:N-1
        x = x + h
        s = s + 2*f(x)
    end
    s = s + f(b)
    return (b-a)*s/(2*N)
end


function integ_trapezoid( f::Array{Float64,1}, a, b, N )
    if length(f) != N + 1
        Base.AssertionError("Length of f is not equal to N + 1")
    end
    h = (b-a)/(N+1)
    s = f[1] + f[N+1]
    #
    for i = 2:N
        s = s + 2*f[i]
    end
    return h*s
end
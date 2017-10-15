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

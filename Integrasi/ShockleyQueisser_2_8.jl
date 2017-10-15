include("integ_trapezoid.jl")

function denumerator(x::Float64)
    return x^3/(exp(x)-1)
end

function numerator(x::Float64)
    return x^2/(exp(x)-1)
end

using PyPlot
const plt = PyPlot

function do_plot_num_denum(a, b, N)
    x = Array( linspace(a, b, N) )
    y1 = zeros(N)
    y2 = zeros(N)
    for i = 1:N
        y1[i] = denumerator(x[i])
        y2[i] = numerator(x[i])
    end
    #
    plt.clf()
    plt.plot(x, y1, marker="o", label="denum")
    plt.plot(x, y2, marker="o", label="numer")
    plt.legend()
    plt.savefig("plt1.png", dpi=300)
end

const SMALL = eps()
const BIG = 10000.0

function u_x_g(x_g)
    A = integ_trapezoid(numerator, x_g, BIG, 30*BIG )
    B = integ_trapezoid(denumerator, SMALL, BIG, 30*BIG )
    return x_g*A/B
end

function do_plot_u_x_g(a, b, N)
    x = Array( linspace(a, b, N) )
    y = zeros(N)
    for i = 1:N
        y[i] = u_x_g(x[i])
    end
    #
    plt.clf()
    plt.plot(x, y, marker="o", label="u_x_g")
    plt.legend()
    plt.savefig("u_x_g.png", dpi=300)
end

function main()
    #println(numerator(SMALL))
    #println(denumerator(SMALL))
    do_plot_u_x_g(SMALL,6.0,100)
end

main()

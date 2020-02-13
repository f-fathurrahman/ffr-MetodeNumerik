function rumus_biasa(a, b, c)
    D = b^2 - 4*a*c
    x1 = (-b + sqrt(D))/(2*a)
    x2 = (-b - sqrt(D))/(2*a)
    return x1, x2
end

function rumus_tidak_biasa(a, b, c)
    D = b^2 - 4*a*c
    x1 = -2*c/(b + sqrt(D))
    x2 = -2*c/(b - sqrt(D))
    return x1, x2
end

function main()
    a::Float64 = 1.0
    b::Float64 = 3000.001
    c::Float64 = -3.0

    println("rumus_biasa = ", rumus_biasa(a,b,c))
    println("rumus_tidak_biasa = ", rumus_tidak_biasa(a,b,c))
end

main()
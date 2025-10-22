# -*- coding: utf-8 -*-
Pkg.activate("ALGFOROPT", shared=true)

using PGFPlots

p = let
    N = 4
    g = GroupPlot(N,1,groupStyle="horizontal sep=0.25cm, ylabels at=edge left",
                      style="width=5cm, xlabel=\$x\$, ylabel=\$y\$, ytick=\\empty")

    f = x->sin(x)

    s = 0.5
    k = 2.0

    a = 7.0
    b = a + s
    ya = f(a)
    yb = f(b)

    x_arr = collect(range(3π/2-2, stop=3π/2+5.5, length=101))

    for i in 1 : N

        ax = Axis([Plots.Linear(x_arr, f.(x_arr), style="solid, black, mark=none"),
              Plots.Scatter([a,b], [f(a), f(b)], style="black, mark=*, mark options={solid, draw=black, fill=black}"),
                ], style="xtick={$a,$b}, xticklabels={\$a\$,\$b\$}, xticklabel style={text height=2ex}")

        push!(g, ax)

        if yb > ya
            a, b = b, a
            ya, yb = yb, ya
            s = -s
        end

        c, yc = b + s, f(b + s)
        if yc > yb
            a, b = a, c
        else
            a, ya, b, yb = b, yb, c, yc
            s *= k
        end
    end

    g
end

plot(p)

typeof(p)

p.axes

save("TEMP_fig.pdf", p)

B = [
888445.0 887112.0;
887112.0 885781.0
]

b = [1.0; 0.0]

x_exact = [885781.0; -887112.0]

B = BigFloat.(B)
b = BigFloat.(b)
x_exact = BigFloat.(x_exact)

display(B); println()
display(b); println()

x = B\b
display(x); println()

display(B*x - b); println() 

display(B*x_exact - b); println()

display(x - x_exact); println()
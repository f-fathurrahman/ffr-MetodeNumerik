def forward_diff_Oh(f, x, h):
    return ( f(x+h) - f(x) )/h

def forward_diff_Oh2(f, x, h):
    return ( -f(x+2*h) + 4*f(x+h) - 3*f(x) )/(2*h)

def backward_diff_Oh(f, x, h):
    return ( f(x) - f(x-h) )/h

def backward_diff_Oh2(f, x, h):
    return ( 3*f(x) - 4*f(x-h) + f(x-2*h) )/(2*h)

def centered_diff_Oh2(f, x, h):
    return ( f(x+h) - f(x-h) )/(2*h)

def centered_diff_Oh4(f, x, h):
    return ( -f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h) )/(12*h)

def f(x):
    return -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.2

x = 0.5
h = 0.25
true_val = -0.9125

print("Using h = ", h)

print()
print("Forward diff")
#
df = forward_diff_Oh(f, x, h)
print("Oh  = %18.10f %18.10e" % (df, abs(df-true_val)) )
#
df = forward_diff_Oh2(f, x, h)
print("Oh2 = %18.10f %18.10e" % (df, abs(df-true_val)) )


print()
print("Backward diff")
#
df = backward_diff_Oh(f, x, h)
print("Oh  = %18.10f %18.10e" % (df, abs(df-true_val)) )
#
df = backward_diff_Oh2(f, x, h)
print("Oh2 = %18.10f %18.10e" % (df, abs(df-true_val)) )

print()
print("Centered diff")
#
df = centered_diff_Oh2(f, x, h)
print("Oh2 = %18.10f %18.10e" % (df, abs(df-true_val)) )
#
df = centered_diff_Oh4(f, x, h)
print("Oh4 = %18.10f %18.10e" % (df, abs(df-true_val)) )
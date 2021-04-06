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
h1 = 0.5
h2 = 0.25
true_val = -0.9125

# Using forward_diff_Oh
Dh1 = forward_diff_Oh(f, x, h1)
Dh2 = forward_diff_Oh(f, x, h2)

#Dh1 = centered_diff_Oh2(f, x, h1)
#Dh2 = centered_diff_Oh2(f, x, h2)

# Richardson extrapolation
Dh12 = 4*Dh2/3 - Dh1/3

print("Dh1  = %18.10f %18.10e" % (Dh1, abs(Dh1-true_val)) )
print("Dh2  = %18.10f %18.10e" % (Dh2, abs(Dh2-true_val)) )
print("Dh12 = %18.10f %18.10e" % (Dh12, abs(Dh12-true_val)) )

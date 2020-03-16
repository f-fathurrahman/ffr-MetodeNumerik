
def integ_trapz( fa, fb, a, b ):
    I = (b - a) * (fa + fb) / 2
    return I

x = [0.0, 0.12, 0.22, 0.32, 0.36, 0.40,
     0.44, 0.54, 0.64, 0.70, 0.80]

fx = [0.200000, 1.309729, 1.305241, 1.743393, 2.074903, 2.456000, 
      2.842985, 3.507297, 3.181929, 2.363000, 0.232000]

I_exact = 1.640533 # from the book

Ndata = len(x)
I = 0.0
for i in range(Ndata-1):
    I = I + integ_trapz( fx[i], fx[i+1], x[i], x[i+1] )

E_t = (I_exact - I)
ε_t = E_t/I_exact * 100
print("Integral result = %.6f" % I)
print("True error      = %.6f" % E_t)
print("ε_t             = %.1f%%" % ε_t)

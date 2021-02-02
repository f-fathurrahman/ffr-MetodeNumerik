from uncertainties import ufloat
from uncertainties.umath import *

F = ufloat(750.0, 30) # N/m
L = ufloat(9, 0.03) # m
E = ufloat(7.5e9, 5e7) # N/m^2
I = ufloat(0.0005, 0.000005) # m^4

y = F*L*L*L*L/(8*E*I)
print("y = ", y)

# Example 6.4 of Kiusalaas
# using trapezoidal rule to evaluate
# \begin{equation}
# \int_{0}^{\pi} \sqrt{x} \cos{x}\,\mathrm{d}x
# \end{equation}

import math
from integ import *

# the function to be integrated
def f(x):
    return math.sqrt(x)*math.cos(x)

Iold = 0.0
kFound = 0
for k in range(1,31):
    Inew = integ_trapezeoid_recursive( f, 0.0, math.pi, Iold, k )
    diffI = abs(Inew-Iold)
    print('k = %5d, Inew = %18.10f, diffI = %18.10e' % (k, Inew, diffI))
    if k > 1 and diffI < 1.0e-6:
        kFound = k
        break
    Iold = Inew

print('Integral result = %18.10f' % Inew)
print('kFound = %d' % kFound)
print('nPanels = %d' % 2**(kFound-1))
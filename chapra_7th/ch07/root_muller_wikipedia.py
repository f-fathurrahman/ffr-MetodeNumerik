from typing import *
from cmath import sqrt  # Use the complex sqrt as we may generate complex numbers

Num = Union[float, complex]
Func = Callable[[Num], Num]


def div_diff(f: Func, xs: List[Num]):
    """Calculate the divided difference f[x0, x1, ...]."""
    if len(xs) == 2:
        a, b = xs
        return (f(a) - f(b)) / (a - b)
    else:
        return (div_diff(f, xs[1:]) - div_diff(f, xs[0:-1])) / (xs[-1] - xs[0])


def mullers_method(f: Func, xs: (Num, Num, Num), iterations: int) -> float:
    """Return the root calculated using Muller's method."""
    x0, x1, x2 = xs
    for _ in range(iterations):
        w = div_diff(f, (x2, x1)) + div_diff(f, (x2, x0)) - div_diff(f, (x2, x1))
        s_delta = sqrt(w ** 2 - 4 * f(x2) * div_diff(f, (x2, x1, x0)))
        denoms = [w + s_delta, w - s_delta]
        # Take the higher-magnitude denominator
        x3 = x2 - 2 * f(x2) / max(denoms, key=abs)
        # Advance
        x0, x1, x2 = x1, x2, x3
        print("x3 = ", x3)
        if abs(f(x3)) < 1e-10:
            print("Converged")
            break
    return x3


def f_example(x: Num) -> Num:
    """The example function. With a more expensive function, memoization of the last 4 points called may be useful."""
    return x ** 2 - 612

import cmath
def chapra_exercise_6_7(x):
    return cmath.sin(x) + cmath.cos(1 + x**2) - 1

root = mullers_method(f_example, (10, 20, 30), 5)
print("Root: {}".format(root))  # Root: (24.738633317099097+0j)

root = mullers_method(chapra_exercise_6_7, (0.5, 1.0, 1.5), 50)
print("Root: {}".format(root))

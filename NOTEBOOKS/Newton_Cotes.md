# General form

Newton-Cotes formula replace a complicated function of tabulated data with an approximating function that is easy to integrate. Specifically, they approximate a definite integral by the following form:
$$
I = \int_{a}^{b} f(x)\,\mathrm{d}x \approx \int_{a}^{b} f_{n}(x)\,\mathrm{d}x
$$
where $f_n(x)$ is a polynomial of the form:
$$
f_{n}(x) = a_{0} + a_{1}x + \cdots + a_{n-1}x^{n-1} + a_{n}x^{n}
$$

# Trapezoidal rule

$$
I = \int_{a}^{b} f(x)\,\mathrm{d}x \approx \int_{a}^{b} f_{1}(x)\,\mathrm{d}x
$$

$$
I = (b-a)\frac{f(a) + f(b)}{2}
$$

Error of trapezoidal rule:
$$
E_{t} = \frac{1}{12}f''(\xi)(b - a)^3
$$
Multiple application of trapezoidal rule:
$$
I = \frac{h}{2} \left[ f(x_0) + 2\sum_{i=1}^{N-1} f(x_i) + f(x_N)
\right]
$$

$$
h = \frac{b-a}{N}
$$

There are $N + 1$ points.
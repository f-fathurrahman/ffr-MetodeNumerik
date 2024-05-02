Time-dependent Schroedinger equation:
$$
\imath \frac{\partial \psi(x,t)}{\partial t} =
\left[ -\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi(x,t)
$$
Spatial mesh:
$$
x_{j} = (j - 1)\Delta x \quad j = 1,2,\ldots,N_{x}
$$
Using central-difference representation of Hamiltonian
$$
\psi^{n+1}_{j} + \imath \frac{\Delta t}{2} \left[
-\frac{\psi_{j+1}^{n+1} - 2\psi_{j}^{n+1} + \psi_{j-1}^{n+1}}{4\Delta x^2} +
V_{j}\psi_{j}^{n+1} \right] =
\psi_{j}^{n} - \imath \frac{\Delta t}{2} \left[
-\frac{\psi_{j+1}^{n} - 2\psi_{j}^{n} + \psi_{j-1}^{n}}{4\Delta x^2} +
V_{j}\psi_{j}^{n} \right]
$$

$$
\lambda = \frac{\Delta t}{4\Delta x}
$$

$$
-\lambda \psi_{j+1}^{n+1} +
\left( 2\lambda + \frac{\Delta t}{2}V_{j} - \imath \right)\psi_{j}^{n+1} -
\lambda \psi_{j-1}^{n+1} =
\lambda \psi_{j+1}^{n} -
\left( 2\lambda + \frac{\Delta t}{2}V_{j} + \imath \right)\psi_{j}^{n} +
\lambda \psi_{j-1}^{n}
$$

$$
W_{j} = 2\lambda + \frac{\Delta t}{2} V_{j}
$$

$$
-\lambda \psi_{j+1}^{n+1} +
\left( W_{j} - \imath \right)\psi_{j}^{n+1} -
\lambda \psi_{j-1}^{n+1} =
\lambda \psi_{j+1}^{n} -
\left( W_{j} + \imath \right)\psi_{j}^{n} +
\lambda \psi_{j-1}^{n}
$$

$$
\psi_{1}^{n+1} = \psi_{N_x}^{n+1} = 0
$$

s

s

1

2

3

4

s

s

s

s

s

s

s


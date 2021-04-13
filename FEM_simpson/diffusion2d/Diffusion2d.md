Heat equation in 2d:
$$
\rho c \frac{\partial T}{\partial x} = -\nabla \cdot \mathbf{q} + A =
\begin{bmatrix}
\dfrac{\partial}{\partial x} & \dfrac{\partial}{\partial y}
\end{bmatrix}
\begin{bmatrix} q_x \\ q_y \end{bmatrix} + A
$$
with constitutive relation (Fourier's law):
$$
\begin{bmatrix} q_x \\ q_y \end{bmatrix} =
-\begin{bmatrix}
k & 0 \\ 0 & k
\end{bmatrix}
\begin{bmatrix}
\dfrac{\partial}{\partial x} \\
\dfrac{\partial}{\partial y} \\
\end{bmatrix} T = -\mathbb{K} \nabla T
$$
to yield
$$
\frac{\partial T}{\partial t} = \nabla (\mathbb{K} \nabla T) + \frac{A}{\rho c}
$$
For constant thermal diffusivity:
$$
\frac{\partial T}{\partial t} = \kappa \left(
\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2}
\right)
$$
Initial condition:
$$
T(x,y,t=0) = 1, \ \ \ \forall x \in [0,Lx], \ \ \ \forall y \in [0,L_y]
$$
Boundary condition:
$$
\begin{align}
T(x=0,y,t) &= 0 \\
T(x=L_x,y,t) &= 0 \\
T(x,y=0,y,t) &= 0 \\
T(x,y=L_y,t) &= 0
\end{align}
$$
Shape functions, defined on local coordinate $[-1,1]\times[-1,1]$ :
$$
N_1(\xi,\eta) = \frac{1}{4} (1 - \xi) (1 - \eta) \\
N_2(\xi,\eta) = \frac{1}{4} (1 - \xi) (1 + \eta) \\
N_3(\xi,\eta) = \frac{1}{4} (1 + \xi) (1 + \eta) \\
N_4(\xi,\eta) = \frac{1}{4} (1 + \xi) (1 - \eta) \\
$$

Approximation of continuous variable $T(x,y)$ as
$$
T \approx \begin{bmatrix}
N_1 & N_2 & N_3 & N_4
\end{bmatrix}
\begin{bmatrix}
T_1 \\ T_2 \\ T_3 \\ T_4
\end{bmatrix} = \mathbf{N}^{T}\mathbf{T}
$$
Residual:
$$
\mathbf{R} = \frac{\partial}{\partial t} \mathbf{N}^{T} \mathbf{T} -
\nabla (\mathbb{K}\nabla \mathbf{N}^{T} \mathbf{T} ) - H
$$
Derivatives of the shape functions in physical coordinates:
$$
\nabla \mathbf{N}^{T} = \begin{bmatrix}
\dfrac{\partial}{\partial x} \\ \dfrac{\partial}{\partial y0}
\end{bmatrix}
\begin{bmatrix}
N_1 & N_2 & N_3 & N_4
\end{bmatrix}
$$

## Integration of Element Matrices

$$
\int_{-1}^{1} f(\xi,\eta)\ \mathrm{d}\xi \mathrm{d}\eta \approx
\sum_{i=1}^{nipx}\sum_{j=1}^{nipy} w_{i} w_{j} f(\xi_{i},\eta_{j}) =
\sum_{k=1}^{nip} w_{k} f(\xi_{k},\eta_{k})
$$


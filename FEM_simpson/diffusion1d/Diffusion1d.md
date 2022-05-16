# Persamaan Difusi 1d

Persamaan difusi:
$$
\frac{\partial T}{\partial t} = \kappa \frac{\partial^2 T}{\partial x^2} + H
$$
dengan kondisi awal:
$$
T(x,t=0) = 0, \ \ \ \forall x \in [0,L_x]
$$
dan kondisi batas:
$$
T(x=0,t) = 0 \\
T(x=L_x,t) = 0
$$
$\kappa$ adalah difusivitas termal

$H$ adalah konstata sumber kalor $10^{-6}$ $\mathrm{K}\mathrm{s}^{-1}$.

$L_x = 10000$ m

## Diskritisasi Elemen Hingga

Persamaan terdiskritisasi:
$$
\mathbf{L}\ \mathbf{T}^{n+1} = \mathbf{R}\ \mathbf{T}^{n} + \mathbf{F}
$$

$$
\mathbf{L} = \frac{\mathbf{M}}{\Delta t} + \mathbf{K}
$$

Mass matrix
$$
\mathbf{M} = \begin{bmatrix}
\dfrac{\Delta x}{3} & \dfrac{\Delta x}{6} \\
\dfrac{\Delta x}{6} & \dfrac{\Delta x}{3}
\end{bmatrix}
$$

$$
\mathbf{R} = \frac{\mathbf{M}}{\Delta t}
$$

Stiffness matrix
$$
\mathbf{K} = \kappa\begin{bmatrix}
\dfrac{1}{\Delta x} & -\dfrac{1}{\Delta x} \\
-\dfrac{1}{\Delta x} & \dfrac{1}{\Delta x}
\end{bmatrix}
$$
Load vector
$$
\mathbf{F} = H\begin{bmatrix}
\dfrac{\Delta x}{2} \\
\dfrac{\Delta x}{2}
\end{bmatrix}
$$

$$
\mathbf{T} = \begin{bmatrix}
T_1 \\ T_2
\end{bmatrix}
$$

 Sistem global (menggunakan matriks dan vektor global):
$$
\tilde{\mathbf{L}} \tilde{\mathbf{T}}^{n+1} = \tilde{\mathbf{R}}\ 
\tilde{\mathbf{T}}^{n} + \tilde{\mathbf{F}} = \tilde{\mathbf{b}}
$$

## Mass matrix

$$
\mathbf{M} = \int_{0}^{L} \mathbf{N}^T\mathbf{N}\ \mathrm{d}x =
\int_{0}^{L} \begin{bmatrix}
N_1 N_1 & N_1 N_2 \\
N_2 N_1 & N_2 N_2
\end{bmatrix}\ \mathrm{d}x
$$

## Stiffness matrix

$$
\mathbf{M} = \int_{0}^{L} \kappa
\frac{\partial \mathbf{N}^T}{\partial x}
\frac{\partial \mathbf{N}}{\partial x} \mathrm{d}x =
\int_{0}^{L} \begin{bmatrix}
\dfrac{\partial N_1}{\partial x} \dfrac{N_1}{\partial x} &
\dfrac{\partial N_1}{\partial x} \dfrac{N_2}{\partial x} \\
\dfrac{\partial N_2}{\partial x} \dfrac{N_1}{\partial x} & 
\dfrac{\partial N_2}{\partial x} \dfrac{N_2}{\partial x}
\end{bmatrix}\ \mathrm{d}x
$$

## Load vector

$$
\mathbf{F} = \int_{0}^{L} H \mathbf{N}^{T}\ \mathrm{d}x =
\int_{0}^{L} H \begin{bmatrix} N_1 \\ N_2 \end{bmatrix} \ \mathrm{d}x
$$

Linear shape functions:
$$
N_{1} = 1 - \frac{x}{L} \\
N_{2} = \frac{x}{L}
$$
Using vector notation:
$$
\mathbf{N} = \begin{bmatrix}
N_1 & N_2
\end{bmatrix}
$$
Shape function derivatives:
$$
\frac{\partial \mathbf{N}}{\partial x} = \begin{bmatrix}
\dfrac{\partial N_1}{\partial x} & \dfrac{\partial N_2}{\partial x}
\end{bmatrix}
$$

## Gauss-Legendre quadrature

Tinjau integral:
$$
\int_{-1}^{1} f(\xi)\ \mathrm{d}\xi
$$
Aturan integral Gauss-Legendre:
$$
\int_{-1}^{1} f(\xi)\ \mathrm{d}\xi \approx
\sum_{k=1}^{n} w_{k} f(\xi_{k})
$$

## Koordinat Lokal

Pemetaan dari $[0,L]$ ke $[-1,1]$:
$$
x = \frac{L}{2}(\xi + 1)
$$
Dengan mensubstitusikan ini ke definisi dari shape function, kita memperoleh definisi shape function pada koordinat lokal:
$$
N_1(\xi) = \frac{1 - \xi}{2} \\
N_2(\xi) = \frac{1 + \xi}{2}
$$
Fungsi pada koordinat lokal dapat digunakan untuk memetakan koordinat dari koordinat lokal ke koordinat ruang fisis. Misalkan jika $x_1$ dan $x_2$ adalah koordinat nodal pada ruang fisis dan $N_1$ dan $N_2$ dievaluasi pada posisi tertentu $\xi$ pada domain $[-1,1]$, maka dengan menggunakan
$$
x = N_1(\xi) x_1 + N_2(\xi) x_{2} = \mathbf{N}\mathbf{x}
$$
memberikan koordinat $x$ pada ruang fisis yang berkorespondensi dengan $\xi$. (dot product?, karena $\mathbf{N}$ adalah vektor baris, $\mathbf{x}$ adalah vektor kolom).

Persamaan ini memiliki bentuk yang sama dengan persamaan yang digunakan untuk mengaproksimasi solusi FEM. Elemen yang memiliki sifat ini disebut sebagai isoparametrik.

Turunan dari shape function dapat didefinisikan pada koordinat lokal.
$$
\frac{\partial N_i}{\partial x} =
\frac{\partial N_i}{\partial \xi}
\frac{\partial \xi}{\partial x}
$$
Dengan menggunakan pemetaan dari $x$ ke $\xi$ diperoleh:
$$
\frac{\partial x}{\partial \xi} = \frac{L}{2}
$$
 Sedangkan:
$$
\frac{\partial N_1}{\partial \xi} = -\frac{1}{2} \\
\frac{\partial N_2}{\partial \xi} = \frac{1}{2}
$$

## Evaluasi integral

Misalkan kita ingin menghitung integral:
$$
\int_{0}^{L} N_{1}(x) N_{1}(x)\ \mathrm{d}x
$$
Pertama, kita mentransformasi integral in ke domain $[-1,1]$ sehingga dapat dievaluasi dengan kuadratur Gauss-Legendre, yang diperoleh dengan mengganti $N_1(x)$ dengan $N_{2}(\xi)$ mensubstitusikan $\mathrm{d}x = (L/2)\mathrm{d}\xi$, sehingga diperoleh
$$
\int_{0}^{L} N_{1}(x) N_{1}(x)\ \mathrm{d}x = \frac{L}{2}
\int_{-1}^{1} N_{1}(\xi) N_{1}(\xi)\ \mathrm{d}\xi
$$
Aplikasi kuadratur Gauss-Legendre:
$$
\frac{L}{2}
\int_{-1}^{1} N_{1}(\xi) N_{1}(\xi)\ \mathrm{d}\xi \approx
\frac{L}{2} \sum_{k=1}^{n} w_{k} N_{1}(\xi_{k}) N_{1}(\xi_{k})
$$
Tinjau integral:
$$
\int_{0}^{L} \kappa \frac{\partial N_2}{\partial x}
\frac{\partial N_2}{\partial x} \ \mathrm{d}x =
\frac{2}{L} \int_{-1}^{1} \kappa \frac{\partial N_2}{\partial \xi}
\frac{\partial N_2}{\partial \xi} \ \mathrm{d}\xi
$$
yang dapat diaproksimasi dengan kuadratur Gauss-Legendre:
$$
\frac{2}{L} \int_{-1}^{1} \kappa \frac{\partial N_2}{\partial \xi}
\frac{\partial N_2}{\partial \xi} \ \mathrm{d}\xi \approx
\frac{2}{L} \sum_{k=1}^{n} \kappa \ w_{k}
\frac{\partial N_2}{\partial \xi_{k}}
\frac{\partial N_2}{\partial \xi_{k}}
$$

## Variable material properties

Tinjau integral:
$$
\int_{0}^{L} H(x) N_{2}(x)\ \mathrm{d}x = \frac{L}{2}\int_{-1}^{1} \approx
\frac{L}{2}\sum_{k=1}^{n} w_{k} H(\xi_{k}) N_{2}(\xi_{k})
$$
dengan mengasumsikan nilai $H_1 = H(\xi_1)$ dan $H_2 = H(\xi_2)$ diketahui, kita dapat menggunaan interpolasi
$$
H(\xi_{k}) = H_1 N_{1}(\xi_{k}) + H_2 N_{2}(\xi_{k})
$$


## Analytic Solution (from Jaeger)

The slab with heat produced within it.
$$
\frac{\partial^2 v}{\partial x^2} - \frac{1}{\kappa}\frac{\partial v}{\partial t} = -\frac{A}{K}
$$
$A$: heat production per unit time per unit volume

$A = A_0$ constant for $t>0$.
$$
v = \frac{A_0 l^2}{2K}\left(
1 - \frac{x^2}{l^2} - \frac{32}{\pi^3} \sum_{n=0}^{\infty}
\frac{(-1)^n}{(2n + 1)^3} \cos \frac{(2n+1)\pi x}{2l}
\exp\left[ -\kappa(2n + 1)^2 \pi^2 t / 4l^2 \right]
\right)
$$
s

s

s

s

s

s

s

s

s

s

s

s

s


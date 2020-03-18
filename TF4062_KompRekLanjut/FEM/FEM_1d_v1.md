# Langkah-langkah Analisis Elemen Hingga

1. Diskritisasi (atau representasi) dari suatu domain menjadi suatu kumpulan elemen hingga.

2. Penurunan persamaan elemen untuk semua elemen pada mesh.

   a. Formulasi variasional dari suatu persamaan diferensial untuk elemen tipikal.

   b. Asumsikan variabel dependen $u$ memiliki bentuk:
   $$
   u = \sum_{i}^{n} u_{i} \psi_{i}
   $$
   dan substitusi ke formulasi variasional untuk mendapat persamaan elemen dalam bentuk:
   $$
   [K^{e}]\{u^{e}\} = \{F^{e}\}
   $$
   
3. TODO

## Model masalah syarat batas

Ingin dicari suatu fungsi $u(x)$ yang memenuhi persamaan diferensial:
$$
-\frac{\mathrm{d}}{\mathrm{d}x}\left(a\frac{\mathrm{du}}{\mathrm{d}x}\right) +
cu - f = 0
$$
untuk domain $0 < x < L$, dengan syarat batas:
$$
u(0) = u_0,\,\,\,\left.\left(a\frac{\mathrm{d}u}{\mathrm{d}x}\right)\right|_{x=L} = Q_{0}
$$
di mana $a = a(x)$, $c = c(x)$, $f = f(x)$, dan $u_0$ dan $Q_{0}$ adalah data (kuantitas yang diketahui dari masalah).

Aproksimasi polinomial dari suatu elemen hingga $\Omega$ diasumsikan dalam bentuk
$$
u^{e}_{h} = \sum_{j=1}^{n} u_{j}^{e} \psi_{j}^{e}(x)
$$
di mana $u_{j}^{e}$ adalah nilai dari solusi $u(x)$ pada node dari elemen hingga $\Omega_{e}$ dan $\psi_{j}^{e}$ adalah fungsi aproksimasi pada elemen hingga.

Bentuk weighted-integral:
$$
0 = \int_{x_a}^{x_b} w\left[
-\frac{\mathrm{d}}{\mathrm{d}x}\left(a \frac{\mathrm{d}u}{\mathrm{d}x} \right) +
cu - f
\right]\,\mathrm{d}x
$$
subdomain $\Omega_{e}=(x_a,x_b)$ dan $w(x)$ adalah fungsi bobot.

Gunakan integrasi-parsial (*integration by parts*):
$$
0 = \int_{x_a}^{x_b} \left(
a \frac{\mathrm{d}w}{\mathrm{d}x} \frac{\mathrm{d}u}{\mathrm{d}x} +
cwu - wf \right)\,\mathrm{d}x -
\left[ w a \frac{\mathrm{d}u}{\mathrm{d}x} \right]_{x_a}^{x_b}
$$
Boundary term:
$$
\left[ w a \frac{\mathrm{d}u}{\mathrm{d}x} \right]_{x_a}^{x_b}
$$
Secondary variable:

Primary variable:

Kondisi:
$$
u_{h}^{e}(x_a) = u_{1}^{e},\,\,\, u_{h}^{e}(x_b) = u_{2}^{e} \\
\left(-a \frac{\mathrm{d}u}{\mathrm{d}x}\right)_{x=x_{a}} = Q_{1}^{e},\,\,\,
\left(a \frac{\mathrm{d}u}{\mathrm{d}x}\right)_{x=x_{b}} = Q_{2}^{e}
$$
Notasi:
$$
0 = \int_{x_a}^{x_b} \left(
a \frac{\mathrm{d}w}{\mathrm{d}x} \frac{\mathrm{d}u}{\mathrm{d}x} +
cwu - wf \right)\,\mathrm{d}x -
w(x_a) Q_1 - w(x_b) Q_2
$$
Blinear form and bilinear form:
$$
B^{e}(w,u) = \int_{x_a}^{x_b}\left(
a \frac{\mathrm{d}w}{\mathrm{d}x} \frac{\mathrm{d}u}{\mathrm{d}x} + cwu
\right)\,\mathrm{d}x \\
l^{e}(w) = \int_{x_a}^{x_b} w\,f\,\mathrm{d}x + w(x_{a})Q_{1} + w(x_{b})Q_{2}
$$
Weak form:
$$
B^{e}(w,u) = l^{e}(w)
$$
Quadratic functional:
$$
\begin{align}
I^{e}(u) & = \frac{1}{2}B^{e}(u,u) - l^{e}(u) \\
& = \frac{1}{2}\int_{x_a}^{x_b}\left[
a\left(\frac{\mathrm{d}u}{\mathrm{d}x}\right)^2 + cu^2
\right]\,\mathrm{d}x - \int_{x_a}^{x_b} u\,f\,\mathrm{d}x - u(x_a)Q_1 - u(x_b)Q_2
\end{align}
$$

## Aproksimasi solusi

Minimum polinomial orde 1:
$$
u_{h}^{e}(x) = c_1^{e} + c_2^{e}x
$$


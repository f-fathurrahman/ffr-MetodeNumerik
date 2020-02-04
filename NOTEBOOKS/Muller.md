Ingat bahwa metode secant mendapatkan estimasi akar dengan cara melakukan proyeksi garis lurus (persamaan linear) pada sumbu $x$ yang melalui dua nilai titik. Metode Muller menggunakan pendekatan yang sama, namun dengan memproyeksikan sebuah parabola (persamaan kuadrat) yang melalui tiga titik. Misalkan parabola tersebut memiliki persamaan sebagai berikut:
$$
f(x) = a(x - x_2)^2 + b(x - x_2) + c
$$
Kita ingin persamaan ini melewati tiga titik: $[x_0, f(x_0)]$, $[x_1, f(x_1)]$ dan $[x_2, f(x_2)]$ atau:
$$
f(x_0) = a(x_0 - x_2)^2 + b(x_0 - x_2) + c \\
f(x_1) = a(x_1 - x_2)^2 + b(x_1 - x_2) + c \\
f(x_2) = a(x_2 - x_2)^2 + b(x_2 - x_2) + c
$$
Koefisien $c$ tidak lain adalah nilai fungsi yang dievaluasi pada tebakan ketiga yaitu $f(x_2)$. Substitusi pada dua persamaan sebelumnya sehingga diperoleh:
$$
f(x_0) - f(x_2) = a(x_0 - x_2)^2 + b(x_0 - x_2) \\
f(x_1) - f(x_2) = a(x_1 - x_2)^2 + b(x_1 - x_2)
$$
Dengan definisi-definisi berikut:
$$
h_0 = x_1 - x_0 \\
h_1 = x_2 - x_1 \\
\delta_0 = \frac{f(x_1) - f(x_0)}{x_1 - x_0} \\
\delta_1 = \frac{f(x_2) - f(x_1)}{x_2 - x_1}
$$
diperoleh:
$$
\begin{align}
(h_0 + h_1)b - (h_0 + h_1)^2 a & = h_0 \delta_0 + h_1 \delta_1 \\
h_1 b - h_1^2 a & = h_{1} \delta_{1}
\end{align}
$$
yang dapat diselesaikan untuk mendapatkan nilai $a$ dan $b$:
$$
\begin{align}
a & = \frac{\delta_1 - \delta_0}{h_1 + h_0} \\
b & = ah_1 + \delta_1 \\
c & = f(x_2)
\end{align}
$$
Setelah diperoleh koefisien-koefisien untuk persamaan kuadrat kita dapat memperoleh akar-akar dari persamaan kuadrat:
$$
x_3 = x_2 + \frac{2c}{b \pm \sqrt{b^2 - 4ac}}
$$
Estimasi error:
$$
\epsilon_a = \left| \frac{x_3 -x_2}{x_3} \right|
$$
Persamaan ini akan menghasilkan sepasang akar. Pada Metode Muller, tanda dari akar dipilih agar sama dengan tanda dari $b$. Pilihan ini akan menghasilkan penyebut yang paling besar sehingga estimasi akar adalah yang paling dekat ke $x_2$. Ketika $x_3$ telah ditentukan proses akan diulangi sehingga perlu ada nilai yang dibuang.

- Jika hanya akar real yang ingin ditemukan, maka dua titik awal yang dipilih adalah yang paling dekat dengan estimasi akar baru $x_3$.
- Jika baik akar real dan kompleks yang ingin ditentukan, maka pendekatan sekuensial yang digunakan, yaitu: $x_1$, $x_2$, dan $x_3$ menggantikan $x_0$, $x_{1}$, dan $x_2$.
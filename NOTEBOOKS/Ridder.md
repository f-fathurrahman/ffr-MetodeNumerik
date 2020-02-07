# Metode Ridder

Asumsikan bahwa akar berada dalam selang $(x_1, x_2)$. Pertama hitung $f_3 = f(x_3)$ di mana $x_3$ adalah titik tengah dari selang, $x_3 = (x_1 + x_2)/2$. Kemudian tinjau suatu fungsi:
$$
g(x) = f(x)e^{(x - x_1)Q}
$$
dengan konstanta $Q$ ditentukan dengan mengharuskan titik-titik $(x_1, g_1)$, $(x_2, g_2)$, dan $(x_3, g_3)$, dengan notasi $g_{i}=g(x_i)$ berada pada garis lurus. Nilai tebakan dari akar diperoleh dengan menggunakan interpolasi linear dari $g(x)$. Dengan menggunakan definisi dari $g(x)$ kita mendapatkan:
$$
g_{1} = f_{1}, \,\,\, g_2 = f_{2}e^{2hQ}, \,\,\, g_{3} = f_{3}e^{hQ}
$$
di mana $h = (x_2 - x_2)/2$. Dengan syarat bawah tiga titik berada pada garis lurus maka $g_3 = (g_1 + g_2)/2$ atau:
$$
f_3 e^{hQ} = \frac{1}{2}\left( f_1 + f_2 e^{2hQ} \right)
$$
yang merupakan persamaan kuadrat dalam $e^{hQ}$. Solusinya adalah:
$$
e^{hQ} = \frac{f_3 \pm \sqrt{f_3^2 - f_1 f_2}}{f_2}
$$
Dengan interpolasi linear berdasarkan titik $(x_1, g_1)$ dan $(x_3, g_3)$ diperoleh nilai tebakan akar:
$$
x_4 = x_3 - g_3\frac{x_3 - x_1}{g_3 - g_1} =
x_3 - f_3 e^{hQ}\frac{x_3 - x_1}{f_3 e^{hQ} - g_1}
$$
Dengan substitusi $e^{hQ}$ diperoleh:
$$
x_4 = x_3 \pm (x_3 - x_1)\frac{f_3}{\sqrt{f_3^2 - f_1 f_2}}
$$
Tanda $(+)$ yang dipilih jika $f_1 - f_2 > 0$ dan tanda $(-)$ bila sebaliknya.
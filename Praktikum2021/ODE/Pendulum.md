Gerak suatu pendulum dapat dinyatakan dengan persamaan diferensial:
$$
\theta''(t) = -\frac{g}{L}\sin(\theta(t)) - k\theta'(t)
$$
dengan syarat awal
$$
\theta(0) = \theta_0
$$
dan
$$
\theta'(0) = 0
$$

Dalam persamaan tersebut $$\theta$$ menyatakan simpangan pendulum,
$$\theta_0$$ menyatakan simpangan awal pendulum,
$$g$$ menyatakan percepatan gravitasi, $L$ menyatakan panjang benang pendulum,
dan $$k\theta'$$ menyatakan suku redaman (gesekan) yang berbanding lurus
dengan kecepatan $$\theta'$$ ($$k$$ diasumsikan merupakan bilangan positif).

Cari solusi $$\theta(t)$$ untuk kasus $$k=0$$ untuk simpangan awal $$\theta_0 = 0.1, \pi/5$$
dan 3.0.

Masih untuk kasus $$k=0$$, tentukan periode osilasi $$T$$ sebagai
fungsi dari simpangan awal $$u_0$$.

Carilah solusi $$\theta(t)$$ untuk kasus $k = 0.1, 0.3, 0.5$ dan 0.8 untuk
simpangan awal yang sama.

Buatlah animasi/video sederhana dengan menggunakan Matplotlib untuk visualisasi gerakan pendulum
yang dijelaskan oleh persamaan diferensial di atas.
Anda dapat menggunakan contoh berikut ini untuk membuat animasi/video sederhana dengan
Matplotlib
http://kuliah2013.tf.itb.ac.id/mod/resource/view.php?id=5689


Gunakan metode Runge-Kutta orde 4.
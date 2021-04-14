Sistem persamaan:
$$
\frac{\partial A}{\partial t} = \nabla^2 A + \gamma(a - A - A^2B) \\
\frac{\partial B}{\partial t} = d\nabla^2 B + \gamma(b - A^2 B)
$$
Sistem persamaan ini digunakan untuk memodelkan interaksi antara dua zat kimia $A$ dan $B$. $B$ menyatakan konsentrasi substrat (inhibitor) yang dikonsumsi pada suatu reaksi oleh suatu aktivator dengan konsentrasi $A$.

Domain: 
$$
A(x,y,t=0) = a + b + r \\
B(x,y,t=0) = \frac{b}{(a + b)^2} + r
$$

$$
\mathbf{M}\frac{\partial}{\partial t} \mathbf{A} + \gamma \mathbf{M}\mathbf{A} +
\mathbf{K}_{A}\mathbf{A} = \mathbf{F}_{A} \\
\mathbf{M}\frac{\partial}{\partial t} \mathbf{B} + \mathbf{K}_{B}\mathbf{B} = \mathbf{F}
$$

s
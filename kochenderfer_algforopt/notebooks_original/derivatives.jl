# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.5.1
#     language: julia
#     name: julia-1.5
# ---

# # Derivatives
#
# This notebook was automatically generated from the Algorithms for Optimization source code. Each cell generates a figure from the original text. While this code is not optimized for use in lectures, we provide it here to be adapted for such projects. We hope you find it useful.

include("support_code.jl");

# +
p = let

f = x->sin(x)

x₀ = 0.5
ymin = -1.1
lo = -2
hi = 3

function get_axis(h)
a, b = x₀, x₀ + h
fa, fb = f(a), f(b)
m = (fb - fa) / (b-a)

p = Plots.Plot[]
push!(p, Plots.Linear(f, (lo,hi), style="thick, black"))
push!(p, Plots.Linear([a,a],[fa, ymin], style="solid, black, mark=none"))
push!(p, Plots.Linear([b,b],[fb, ymin], style="solid, black, mark=none"))
push!(p, Plots.Linear([lo,hi],[ymin, ymin], style="solid, black, mark=none"))
push!(p, Plots.Linear([lo, hi],[fa + (lo-a)*m, fa + (hi-a)*m], style="ultra thick, solid, pastelBlue, mark=none"))
push!(p, Plots.Scatter([a, b],[fa, fb], style="solid, mark=*, mark options={draw=black, fill=black}"))

Axis(p, style="axis lines=left, axis line style = {transparent}, xtick={$(x₀ + h/2)}, xticklabels={\$h\$}, ytick=\\empty, xtick style={draw=none}, ytick style={draw=none}") #$
end

g = GroupPlot(3, 1, groupStyle="horizontal sep=1cm", style="ymin=-1.1, ymax=1.1, width=6cm")
push!(g, get_axis(2.0))
push!(g, get_axis(1.0))
push!(g, get_axis(0.5))
g
end

plot(p)

# +
% 	plot(Axis(
% 		[Plots.Contour(x->exp(0-x[1]^2-x[2]^2)*(x[1]-2*x[2]), (-2,2), (-2,2), levels=[-0.1,-0.3, -0.5]),
% 		 Plots.Scatter([-0.75], [1.0], style="black, mark=*, mark options={draw=black}"),
% 		 Plots.Command("\\draw [->, ultra thick] (axis cs:-0.75,1) -- (axis cs:-0.25,0.5);"),
% 		 Plots.Command("\\node at (axis cs:-0.75,1) [anchor=north] {\$\\vect x\$};"),
% 		 Plots.Command("\\node at (axis cs:-0.5,0.75) [anchor=south west] {\$\\vect s\$};"),
% 		],
% 	          style="width=9cm, axis lines=left, xtick=\\empty, ytick=\\empty, xmin=-2, xmax=1.75, ymin=-0.5, ymax=2, view={0}{90}, contour/labels=false, contour/draw color={black}"))
% \end{jlcode}
% \begin{marginfigure}[-9cm]
% 	\centering
% 	\plot{fig_derivatives_2d_derivative}
% 	\caption{The two-dimensional directional derivative}
% \end{marginfigure}

\begin{example}
	We wish to compute the directional derivative of $f(\vect x) = x_1 x_2$ at $\vect x = [1,0]$ in the direction $\vect s = [-1, -1]$:
	\begin{align*}
		\nabla f(\vect x) & = \begin{bmatrix}\frac{\partial f}{\partial x_1}, & \frac{\partial f}{\partial x_2}\end{bmatrix} = \brock*{x_2, x_1} \\
		\nabla_{\vect s}f(\vect x) & = \nabla f(\vect x)^\transpose \vect s = \begin{bmatrix}0 & 1\end{bmatrix} \begin{bmatrix}-1 \\ -1\end{bmatrix} = -1
	\end{align*}
	We can also compute the directional derivative as follows:
	\begin{align*}
		g(\alpha) &= f(\vect x + \alpha \vect s) = (1 - \alpha)(-\alpha) = \alpha^2 - \alpha \\
		g'(\alpha) &= 2\alpha - 1 \\
		g'(0) &= -1
	\end{align*}
	\caption{Computing a directional derivative.\label{ex:directional-deriv}}
\end{example}
\begin{marginfigure}[-8cm]
	\centering
	\begin{tikzpicture}
	\begin{axis}[
	width=\marginparwidth,
	xlabel={$\alpha$},
	ylabel={$f(\alpha)$},
	axis lines=middle, xtick={1}, ytick=\empty,
	]
	\addplot[mark=none, samples=100, domain=-0.25:1.25]{(1-x)*x};
	\addplot[mark=none, samples=2, domain=-0.25:0.5]{x};
	\end{axis}
	\end{tikzpicture}
\end{marginfigure}



\begin{ignore}
	For multi-input, multi-output functions\sidenote{$\vect f: \mathbb{R}^n \rightarrow \mathbb{R}^m$}, we can compute the \vocab{Jacobian}:

	\begin{equation}
	\nabla \vect f(\vect x) = \begin{bmatrix}
	\frac{\partial f_1(\vect x)}{\partial x_1} & \frac{\partial f_1(\vect x)}{\partial x_2} & \cdots & \frac{\partial f_1(\vect x)}{\partial x_n} \\
	& \vdots & & \\
	\frac{\partial f_m(\vect x)}{\partial x_1} & \frac{\partial f_m(\vect x)}{\partial x_2} & \cdots & \frac{\partial f_m(\vect x)}{\partial x_n} \\
	\end{bmatrix}
	\end{equation}

	The Jacobian forms the best linear approximation for $\vect f$ at $\vect x$, because it contains all of the gradient components.
	\tim{I don't think we use the Jacobian in the text. Should we leave it out?} \mykel{Let's comment it out. Let's try to make the book as simple as possible (but no simpler). If a topic is not connected with anything, we should strongly consider dropping it. Also, the Jacobian gets weird when we have the convention where the gradient is a column vector.}
\end{ignore}

\section{Numerical Differentiation}

The process of estimating derivatives numerically is referred to as \vocab{numerical differentiation}.
Estimates can be derived in different ways from function evaluations. This section discusses finite difference methods and the complex step method.\sidenote{For a more comprehensive treatment of the topics discussed in the remainder of this chapter, see \fullcite{Griewank2008}.}


% The gradient $\nabla f$ or the derivative along a search direction ${\frac{\partial}{\partial \alpha} f(\vect x + \alpha \vect d)}$ are often not readily avaiable for direct use, but they are critical to many optimization methods.
% Given a method for derivative estimation one also automatically has a method for gradient estimation - one need merely estimate the partial derivative along each component direction.

\subsection{Finite Difference Methods}

As the name implies, \vocab{finite difference methods} compute the difference between two values that differ by a finite step size.
They approximate the derivative definitions in \cref{eq:derivative_defintions} using small differences:
\begin{fullwidth}
\begin{equation}
	f'(x) \approx \underbrace{\frac{f(x + h) - f(x)}{h}}_\text{forward difference} \approx \underbrace{\frac{f(x + h/2) - f(x - h/2)}{h}}_\text{central difference} \approx \underbrace{\frac{f(x) - f(x - h)}{h}}_\text{backward difference}
\end{equation}
\end{fullwidth}

Mathematically, the smaller the step size $h$, the better the derivative estimate.
Practically, values for $h$ that are too small can result in numerical cancellation errors.
This effect is shown later in \cref{fig:derivatives_derivative_error}. \Cref{alg:finite-diff} provides implementations for these methods.

\begin{algorithm}
\begin{juliaverbatim}
diff_forward(f, x; h=sqrt(eps(Float64))) = (f(x+h) - f(x))/h
diff_central(f, x; h=cbrt(eps(Float64))) = (f(x+h/2) - f(x-h/2))/h
diff_backward(f, x; h=sqrt(eps(Float64))) = (f(x) - f(x-h))/h
\end{juliaverbatim}
	\caption{
		\label{alg:finite-diff}
		Finite difference methods for estimating the derivative of a function \jlv{f} at \jlv{x} with finite difference \jlv{h}.
		The default step sizes are the square root or cube root of the machine precision for floating point values.
		These step sizes balance machine round-off error with step size error.

		The \jlv{eps} function provides the step size between \jlv{1.0} and the next larger representable floating-point value.
	}
\end{algorithm}


The finite difference methods can be derived using the Taylor expansion.
We will derive the forward difference derivative estimate, beginning with the Taylor expansion of $f$ about $x$:
\begin{equation}
	f(x + h) = f(x) + \frac{f'(x)}{1!} h + \frac{f''(x)}{2!} h^2 + \frac{f'''(x)}{3!} h^3 + \cdots
\end{equation}

We can rearrange and solve for the first derivative:
\begin{align}
		f'(x) h & = f(x+h) - f(x) - \frac{f''(x)}{2!} h^2 - \frac{f'''(x)}{3!} h^3 - \cdots \\
		f'(x) & = \frac{f(x+h) - f(x)}{h} - \frac{f''(x)}{2!} h - \frac{f'''(x)}{3!} h^2 - \cdots \\
		f'(x) & \approx \frac{f(x+h) - f(x)}{h}
\end{align}

The forward difference approximates the true derivative for small $h$ with error dependent on $\frac{f''(x)}{2!} h + \frac{f'''(x)}{3!} h^2 + \cdots$.
This error term is $O(h)$, meaning the forward difference has linear error as $h$ approaches zero.\sidenote[][]{
	Asymptotic notation is covered in \chref{the appendix}{ch:math_concepts}.
}

The central difference method has an error term of $O(h^2)$.\cite{Mathews2004}
We can derive this error term using the Taylor expansion.
The Taylor expansions about $x$ for $f(x+h/2)$ and $f(x-h/2)$ are:
\begin{align}
		f(x+h/2) & = f(x) + f'(x) \frac{h}{2} + \frac{f''(x)}{2!} \paren*{\frac{h}{2}}^2 + \frac{f'''(x)}{3!}\paren*{\frac{h}{2}}^3 + \cdots \\
		f(x-h/2) & = f(x) - f'(x) \frac{h}{2} + \frac{f''(x)}{2!} \paren*{\frac{h}{2}}^2 - \frac{f'''(x)}{3!}\paren*{\frac{h}{2}}^3 + \cdots
\end{align}

Subtracting these expansions produces:
\begin{equation}
	f(x+h/2) - f(x-h/2) \approx 2f'(x)\frac{h}{2} + \frac{2}{3!} f'''(x) \paren*{\frac{h}{2}}^3
\end{equation}

We rearrange to obtain:
\begin{equation}
	f'(x) \approx \frac{f(x+h/2) - f(x-h/2)}{h} - \frac{f'''(x)h^2}{24}
\end{equation}
which shows that the approximation has quadratic error.


\subsection{Complex Step Method}

We often run into the problem of needing to choose a step size $h$ small enough to provide a good approximation but not too small so as to lead to numerical subtractive cancellation issues.
The \vocab{complex step method} bypasses the effect of subtractive cancellation by using a single function evaluation. We evaluate the function once after taking a step in the imaginary direction.\sidenote{\fullcite{Martins2003}. Special care must be taken to ensure that the implementation of $f$ properly supports complex numbers as input.}

The Taylor expansion for an imaginary step is:
\begin{equation}
	f(x + i h) = f(x) + i h f'(x) - h^2 \frac{f''(x)}{2!} - i h^3 \frac{f'''(x)}{3!} + \cdots
\end{equation}
Taking only the imaginary component of each side produces a derivative approximation:
\begin{align}
		\Imag(f(x + i h)) & = h f'(x)  - h^3 \frac{f'''(x)}{3!} + \cdots \\
		\Rightarrow f'(x) & = \frac{\Imag (f(x + i h))}{h} + h^2 \frac{f'''(x)}{3!} - \cdots \\
						  & = \frac{\Imag (f(x + i h))}{h} + O(h^2) \text{ as } h \rightarrow 0
\end{align}
An implementation is provided by \cref{alg:complex-step}. The real part approximates $f(x)$ to within $O(h^2)$ as $h \rightarrow 0$:
\begin{align}
		\Real(f(x + i h)) = f(x) - h^2 \frac{f''(x)}{2!} + \hdots \\
		\Rightarrow f(x) = \Real(f(x + i h)) + h^2 \frac{f''(x)}{2!} - \cdots
\end{align}

Thus, we can evaluate both $f(x)$ and $f'(x)$ using a single evaluation of $f$ with complex arguments. \Cref{ex:complex-step} shows the calculations involved for estimating the derivative of a function at a particular point.
\Cref{alg:complex-step} implements the complex step method.
\Cref{fig:derivatives_derivative_error} compares the numerical error of the complex step method to the forward and central difference methods as the step size is varied.


\begin{algorithm}
\begin{juliaverbatim}
diff_complex(f, x; h=1e-20) = imag(f(x + h*im)) / h
\end{juliaverbatim}
\caption{
	The complex step method for estimating the derivative of a function \jlv{f} at \jlv{x} with finite difference \jlv{h}.
	\label{alg:complex-step}
}
\end{algorithm}

# %Use of the complex step method requires that the function be able to accept complex values as input.
# %Existing packages sometimes require modifying the code that implements $f$, which can limit the method's direct application.

\begin{example}
Consider $f(x) = \sin(x^2)$.
The function value at $x = \pi/2$ is approximately \num{0.624266} and the derivative is $\pi \cos(\pi^2/4) \approx -2.45425$.
We can arrive at this using the complex step method:
\begin{juliaconsole}
f = x -> sin(x^2);
v = f(π/2 + 0.001im);
real(v) # f(x)
imag(v)/0.001 # f'(x)
\end{juliaconsole}
\caption{The complex step method for estimating derivatives.\label{ex:complex-step}}
\end{example}

\begin{jlcode}
	diff_forward(f, x; h=sqrt(eps(Float64))) = (f(x+h) - f(x))/h
	diff_central(f, x; h=cbrt(eps(Float64))) = (f(x+h/2) - f(x-h/2))/h
	diff_backward(f, x; h=sqrt(eps(Float64))) = (f(x) - f(x-h))/h
	diff_complex(f, x; h=1e-20) = imag(f(x + h*im)) / h

	p = let
		x = 0.5
		d_true = 1.0

		arr_h = collect(10 .^ range(-17, stop=1, length=101))
		arr_forward = [abs(d_true - diff_forward(sin, x, h=h)/cos(x))
						for h in arr_h]
		arr_central = [abs(d_true - diff_central(sin, x, h=h)/cos(x))
						for h in arr_h]
		arr_complex = [abs(d_true - diff_complex(sin, x, h=h)/cos(x))
						for h in arr_h]

		Axis([
			Plots.Linear(arr_h, arr_complex, style="thick, solid,
				pastelGreen, mark=none", legendentry="complex"),
			Plots.Linear(arr_h, arr_forward, style="thick, solid,
				pastelBlue, mark=none", legendentry="forward"),
			Plots.Linear(arr_h, arr_central, style="thick, solid,
				pastelRed, mark=none", legendentry="central")],
			xmode="log", ymode="log",
			xlabel="step size \$h\$",
			ylabel="absolute relative error",
			width="10cm", height="8cm",
			style="legend pos=outer north east"
		)
	end
	plot(p)

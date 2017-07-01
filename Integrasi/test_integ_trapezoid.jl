include("integ_trapezoid.jl")

function f( t::Array{Float64} )
  return sqrt(9.81*68.1/0.25)*tanh.(sqrt(9.81*0.25/68.1)*t)
end

function f( t::Float64 )
  return sqrt(9.81*68.1/0.25)*tanh(sqrt(9.81*0.25/68.1)*t)
end

using PyPlot
const plt = PyPlot

function do_plot(a, b, N)
  x = Array( linspace(a, b, N) )
  y = f(x)
  #
  plt.clf()
  plt.plot(x, y, marker="o")
  plt.savefig("plt1.png", dpi=300)
end

function main()
  const a = 0.0
  const b = 3.0

  s_old = 0.0
  for i = 1:30
    N = i*20
    s = integ_trapezoid( f, a, b, N )
    @printf("%8d %18.10f %18.10e\n", N, s, abs(s-s_old))
    s_old = s
  end

  println("Program ended normally")
end

main()

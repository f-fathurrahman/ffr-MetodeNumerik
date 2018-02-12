function f = soal_09(Vm)
  q = 1.6022e-19
  kB = 1.3806e-23
  Voc = 0.5
  T = 297
  kBT = kB*T
  RHS = exp(q*Voc/kBT)
  LHS = exp(q*Vm/kBT)*(1 + q*Vm/kBT)
  f = RHS - LHS
endfunction

function eval_LHS(Vm)
  q = 1.6022e-19
  kB = 1.3806e-23
  Voc = 0.5
  T = 297
  kBT = kB*T
  RHS = exp(q*Voc/kBT)  
endfunction

function do_plot()
  Vm_1 = 0.0
  Vm_2 = 1.0
  Npoints = 100
  Vm = linspace(Vm_1, Vm_2, 100)
  f = zeros(1,Npoints)
  for i = 1:Npoints
    f(i) = soal_09(Vm(i))
    printf("%18.10f %18.10f\n", Vm(i), f(i))
  end
endfunction

//do_plot()

exec("bisection.sce", -1)
root = bisection( soal_09, 0.4242424242, 0.4343434343, NiterMax=50 )

exec("ridder.sce", -1)
root = ridder( soal_09, 0.4242424242, 0.4343434343, NiterMax=50 )

if getscilabmode() ~= "STD"
  quit()
end
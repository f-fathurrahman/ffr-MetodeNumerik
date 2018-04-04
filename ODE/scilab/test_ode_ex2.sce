exec("ode_RK4.sce",-1)
exec("ode_ABM.sce",-1)
exec("ode_hamming.sce",-1)

function x = df(t,y)
  x = -y + 1
endfunction

function y = x_analytic(t)
  y = 1 - exp(-t)
endfunction

t0 = 0
tf = 10
y0 = 0
N = 50

KC = 0
for KC = 0:1
  tic, [t1,yR] = ode_RK4(df,[t0 tf],y0,N); tR = toc
  tic, [t1,yA] = ode_ABM(df,[t0 tf],y0,N,KC); tA = toc
  tic, [t1,yH] = ode_hamming(df,[t0 tf],y0,N,KC); tH = toc
  yt = x_analytic(t1)
  //
  NN = length(yt)
  err1 = norm(yt-yR)/NN
  err2 = norm(yt-yA)/NN
  err3 = norm(yt-yH)/NN

  if KC == 1
    printf("With modifier:\n")
  else
    printf("Without modifier:\n")
  end

  printf("RK4     error = %18.10e\n",err1)
  printf("ABM     error = %18.10e\n",err2)
  printf("Hamming error = %18.10e\n",err3)

end

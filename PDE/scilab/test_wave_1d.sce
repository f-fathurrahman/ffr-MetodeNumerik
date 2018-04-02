exec("wave_1d.sce",-1)
exec("to_string.sce",-1)

function y = it0(x)
  y = x.*(1-x)
endfunction

function y = i1t0(x)
  y = 0
endfunction

function y = bx0t(t)
  y = 0
endfunction

function y = bxft(t)
  y = 0
endfunction

a = 1
xf = 1
M = 20
T = 2
N = 50;
[u,x,t] = wave_1d(a,xf,T,it0,i1t0,bx0t,bxft,M,N);

for n = 1:N
  clf()
  plot(x,u(:,n))
  set(gca(),"data_bounds", [0 xf -0.3 0.3])
  xs2png(gcf(), "TEMP_wave_" + to_string(n) + ".png")
end
exec("ode_euler.sce",-1)
exec("ode_euler_PC.sce",-1)
exec("ode_RK4.sce",-1)

// f and y are array with size 2
// y(1) = y
// y(2) = y'
//
// f(1) = y' = y(2)
// f(2) = y'' = -y = -y(1)
// There is no dependence of f to t explicitly.
// However, we put it here because in general there might be dependence
// to t
function f = dy(t,y)
  f(1) =  y(2)
  f(2) = -y(1)
  f = f'
endfunction

tspan = [0 100]
y0 = [0 1]

method = "RK4"

h = 0.05
N = (tspan(2) - tspan(1))/h

if method == "RK4"
  [t,y] = ode_RK4(dy,tspan,y0,N)
elseif method == "euler"
  [t,y] = ode_euler(dy,tspan,y0,N)
elseif method == "euler_PC"
  [t,y] = ode_euler_PC(dy,tspan,y0,N)  
else
  error("method is unknown")
end

clf()
plot( y(:,1), y(:,2), 'b' )
// xmin, ymin, xmax, ymax
xlabel('$y_1$')
ylabel('$y_2$')
if method == "euler"
  square(-15,-15,15,15)
  xs2pdf( gcf(), "images/soal_01_ode_euler_y1_y2.pdf" )
elseif method == "euler_PC"
  square(-1.5,-1.5,1.5,1.5)
  xs2pdf( gcf(), "images/soal_01_ode_euler_PC_y1_y2.pdf" )  
elseif method == "RK4"
  square(-1.5,-1.5,1.5,1.5)
  xs2pdf( gcf(), "images/soal_01_ode_RK4_y1_y2.pdf" )
end

clf()
plot( t, y(:,1), 'b')
xlabel('$t$')
ylabel('$y_1$')
if method == "euler"
  xs2pdf( gcf(), "images/soal_01_ode_euler_t_y1.pdf")
elseif method == "euler_PC"
  set(gca(),"data_bounds",[0,100,-1.2,1.2])
  xs2pdf( gcf(), "images/soal_01_ode_euler_PC_t_y1.pdf")  
elseif method == "RK4"
  set(gca(),"data_bounds",[0,100,-1.2,1.2])  
  xs2pdf( gcf(), "images/soal_01_ode_RK4_t_y1.pdf")
end

if getscilabmode() ~= "STD"
  quit()
end


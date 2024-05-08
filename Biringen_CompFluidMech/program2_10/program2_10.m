%          ***********   PROGRAM 2.10 Nonlinear Convection Equation   *********

%   Solving the nonlinear convection equation using the MacCormack Explicit Method
%    Refer to Biringen and Chow, "Computational Fluid Mechanics by Example," Eqns 2.13.9-2.13.10.
clear all;
umax=10;
xmin=0;
xmax=40;
L=xmax-xmin;
Nx=200;
dx=L/(Nx-1);
C=1.0;
tend=4;
n=1;
t(n)=0;
x=linspace(0,xmax,Nx);
while t(n)<=tend
    if t(n)==0
    for i=1:Nx 
     if x(i)<20
     u(1,i)=umax;
     else
     u(1,i)=0;
      end
      end 
      uhat=u(1,:);
      uin=uhat;
    else
       
    u(n,1)=1;
    for i=2:Nx-1
     fi=1/2*(u(n-1,i))^2;
     fip1=1/2*(u(n-1,i+1))^2;
     uhat(i)=u(n-1,i)-(dt/dx)*(fip1-fi);
      fbi=1/2*(uhat(i))^2;
      fbim1=1/2*(uhat(i-1))^2;
       u(n,i)=1/2*(u(n-1,i)+uhat(i)-(dt/dx)*(fbi-fbim1));
      end
     %Set u = umax Dirichlet boundary condition at x=xmin
     %Boundary condition at umax is u = 0; also Dirichlet
       u(n,Nx)=0;
      u(n,1)=umax;
    end
    A=linspace(1,Nx,1);
    for i =1:Nx
    A(i)=u(n,i)';
    end
    umax=max(A);
    dt=C*dx/umax;
    n=n+1;
    t(n)=t(n-1)+dt;
    
end
%choose time step  to plot, e.g. time step 100.


plot(x,u(100,:),'--');
hold on
plot(x,uin);
title('MacCormack method; C = 1; maximum velocity 10')
axis([0 50 0 12])
xlabel('x')
ylabel('u')
legend(' t=0; --, t=2.01')
%gtext('t = 0')
%gtext('t = 2.01 s')




        
        
                
            
                

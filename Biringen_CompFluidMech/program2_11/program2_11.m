%     ********  PROGRAM 2.11   ****************

% Solving the nonlinear nonlinear convection equation using the Beam & Warming Implicit Method
%Refer to Biringen and Chow "An Introduction to Computational Fluid Mechanics by Eample," Eqns. 
% 2.13.9 - 2.13.19
clear all; close all; clc
umax=10;
XL=0;
XH=40;
NX=101;
L=XH-XL
dx=L/(NX-1);
C=1.0;
tend=4;
n=1;
t(n)=0.;
mu=0.8;
x=linspace (0,XH,NX);
%Initial conditions
for j=1:NX
    %if x(j)<=1
    if x(j)<20
        u(1,j)=umax;
    else
        u(1,j)=0;
    end
end
dt=C*dx/umax
while t(n)<=tend
    LD=-(1/4*dt/dx)*u(n,2:NX-2);
    MD=ones(1,NX-2);
    UD=(1/4*dt/dx)*(u(n,3:NX-1));
    A=diag(LD,-1)+diag(MD,0)+diag(UD,1);
    RHS(1)=u(n,2)+(1/4*dt/dx)*(u(n,1))^2;
    for j=2:NX-3
        artvis=u(n,j+3)-4*u(n,j+2)+6*u(n,j+1)-4*u(n,j)+u(n,j-1);
        RHS(j)=u(n,j+1)-(mu/8)*artvis;
    end
    RHS(j+1)=u(1,NX-1)-(1/4*dt/dx)*(u(n,NX))^2;
    uint=A\RHS';
    %Homogeneous Dirichlet boundary conditions
    u(n+1,1)=u(n,1);
    u(n+1,2:NX-1)=uint';
    u(n+1,NX)=u(n,NX);
    dt=C*dx/max(u(n+1,:));
    n=n+1;
   t(n)=t(n-1)+dt;
end
plot(x,u(1,:));
hold on
plot(x,u(50,:));
plot(x,u(90,:));
%gtext('t = 0')
%gtext('t = 2s')
%gtext('t = 3.6s')
axis([0 45 -2 12])

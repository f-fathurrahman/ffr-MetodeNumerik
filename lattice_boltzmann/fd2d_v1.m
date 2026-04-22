% Finit-Difference 2D, diffusion equation
clear all; close all;
m=51;n=51;
xl=1.0;yl=1.0;
dx=xl/(m-1.0); dy=yl/(n-1.0);
dx2=dx*dx;
dy2=dy*dy;
dt=0.0001;
alpha=0.25;
Twall=1.0;
T=zeros(m,n);To=zeros(m,n);Tm=zeros(m);z=zeros(n,m);
x=zeros(m);y=zeros(n);
x(1)=0.0; y(1)=0.0;
%boundary conditions
for i=1:m-1
x(i+1)=x(i)+dx;
end
for j=1:n-1
y(j+1)=y(j)+dy;
end
%boundry condition, left vertical boundary
for j=1:m
T(1,j)=Twall;
To(1,j)=Twall;
end
nstep=1600;
for kk=1:nstep
for j=2:n-1
for i=2:m-1
xterm=(To(i+1,j)+To(i-1,j))/dx2;
yterm=(To(i,j+1)+To(i,j-1))/dy2;
dd=1./dx2+1./dy2;
T(i,j)=To(i,j)+dt*alpha*(xterm+yterm-2.*To(i,j)*dd);
end
end
%update
for j=2:n-1
for i=2:m-1
To(i,j)=T(i,j);
end
end
%boundary conditions
%bottom, adiabatic
for i=2:m-1
To(i,1)=To(i,2);
end
%right boundary condition T=0
for j=1:n
To(m,j)=0.0;
end
end
for i=1:n
Tm(i)=To(i,(n-1)/2);
end
%matrix rotation for contour plotting
for j=1:n
for i=1:m
z(j,i)=To(i,j);
end
end
figure(1)
plot(x,Tm,'x');
figure (2)
contour(z)
title('Temperature')
xlabel('X')
ylabel('Y')
title('Flux')
xlabel('X')
ylabel('Flux')
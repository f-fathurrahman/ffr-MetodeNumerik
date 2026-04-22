% Chapter 5
% LBM- 1-Finite Difference, flux boundary condition
clear all; close all;
m=101;
dx=1.0;
dx2=dx*dx;
dt=0.5;
qf=100.; tk=20.;
T=zeros(m);To=zeros(m);
x=zeros(m);fluxq=zeros(m);
x(1)=0.0;
T(1)=T(2)+qf*dx/tk; To(1)=T(1);
for i=1:m-1
x(i+1)=x(i)+dx;
end
alpha=0.25;
nstep=800;
for kk=1:nstep
for i=2:m-1
T(i)=To(i)+dt*alpha*(To(i+1)-2.*To(i)+To(i-1))/dx2;
end
%boundary condition
T(1)=T(2)+qf*dx/tk;
%update
for k=1:m-1
To(k)=T(k);
end
To(m)=T(m-1);
end
%Flux:
for k=1:m-1
fluxq(k)=(To(k)-To(k+1)/dx);
end
fluxq(m)=fluxq(m-1);
figure(1)
plot(x,To)
title('Temperature, nstep=400')
xlabel('X')
ylabel('T')
figure(2)
plot(x,fluxq,'x')
title('Flux, time step=400')
xlabel('X')
ylabel('Flux')
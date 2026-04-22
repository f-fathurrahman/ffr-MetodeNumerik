% Finite Difference, 1-D diffusion equation with variable alpha
clear
m=101;
dx=1.0;
dx2=dx*dx;
dt=0.125;
Twall=1.0;
T=zeros(m);To=zeros(m);
x=zeros(m);fluxq=zeros(m);
tk=zeros(m);cpr=zeros(m);
alpha=zeros(m);
x(1)=0.0;
T(1)=Twall; To(1)=Twall;
for i=1:m-1
x(i+1)=x(i)+dx;
end
for i=1:m
tk(i)=20.+30.0/(2.*x(i)+1.);
alpha(i)=tk(i)/100.;
end
nstep=4000;
for kk=1:nstep
for i=2:m-1
tke=0.5*(alpha(i+1)+alpha(i));
tkw=0.5*(alpha(i)+alpha(i-1));
T(i)=To(i)+dt*(tke*(To(i+1)-To(i))-tkw*(To(i)-To(i-1)))/dx2;
end
%update
for k=2:m-1
To(k)=T(k);
end
To(m)=T(m-1);
end
%Flux:
for k=1:m-1
fluxq(k)=tk(k)*(To(k)-To(k+1))/dx;
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
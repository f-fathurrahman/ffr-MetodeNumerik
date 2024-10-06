% PURPOSE: Driver for the one-dimensional (in z-direction)
%  DEM-simulation DEMround1D.m
%  of round particles (finite radius rad). Collisions are 
%  purely elastic: No dissipation has been specified.
% CAVEAT: There may be numerical dissipation, i.e. the
%  total energy is not conserved, because the integrator
%  used (ode113) is not energy-conserving
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

clear all
format compact
n_part=5
% initialize radius and mass
global rad, rad(1:n_part)=0.5;
global m, m(1:n_part)=1;
global E, E=1000; % Young’s modulus
global lmax, lmax=2*n_part+2;
global lmin, lmin=0;
global g, g=-9.81;
% initialize positions and velocities=0
r0=2*[1:n_part];
v0=r0*0;
y0(1:2:2*n_part-1)=r0;
y0(2:2:2*n_part)=v0;
t_end=4
[t,y]=ode113('DEMround1D',[0 t_end],y0);
hold on
for i=1:n_part
  plot(t,y(:,2*i-1),'ko-')
end
axis([0 max(t) lmin-.5 lmax+.5])
return

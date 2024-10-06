function [dy] = linosc_symp_visc_dry(t,y)
% PURPOSE: Returns dy/dt = f(t,y) for the 
%  linear oscillator with no, viscous and dry friction.
% USAGE: To be called from ch1p32symp_visc_dry
% ALGORITHM: 
%  The second order diferential equation
%      ..   .
%    m*x +D*x+k*x=0
%  is transformed into two first-order equation
%    ..  
%    x =(-D*v-k*x)/m
%    . 
%    x =v
%  The input vector y contains the velocities and 
%  positions, y =(v x) which are then used as
%    dy(1) = (-D*y(1)-k*y(2))/m
%    dy(2) = y(1)
% CAVEAT: For large coefficients of dry friction
%  compared to the initial energy, the time step
%  is reduced considerably to obtain velocities close
%  to zero due to the alternation of +-mu as force
%  due to the velocity changes. The resulting reduction
%  of the timestep and high rejection rate for the
%  adaptive algorithms leads to large computation times.
% REVISION HISTORY: 23-May-2014 H.-G. Matuttis
%
global mass   
global k
global D     
global mu     

j=1;
i=2*j-1;
v=y(i);
force=-mu*sign(v)- 2*D*y(i)-k*y(i+1);
dy(i,1)   = force/mass;
dy(i+1,1) = y(i) ;

return
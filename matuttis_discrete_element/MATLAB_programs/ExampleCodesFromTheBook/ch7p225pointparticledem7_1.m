% PURPOSE: Discrete element simulation of point particle
%  bouncing on a floor (at height 0) without
%  dissipation (only elastic forces). The force between 
%  floor and particle is proportional to spring constant
%  k and the penetration into the floor (coordinate).
% USAGE: If the program is used as it is, the height is 
%  plotted at the time where the integrator choose to
%  evaluate the trajectory; if the lower part is commented
%  out, the solution will be interpolated ever at times which
%  are multiples of 0.1
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

%
clear, format compact
global k, k=100;
global g, g=-9.81;
x0=4
v0=1
tspan=[0 10]
[t,y]=ode23('bouncing_ball',tspan,[v0 x0]);
plot(t,y(:,2),'k*')
% uncomment to obtain the continuous trajectories:
%hold on
%tspan=[0:.1:10]
%[t2,y2]=ode23('bouncing_ball',tspan,[v0 x0]);
%plot(t2,y2(:,2),Åfk-Åf)
return
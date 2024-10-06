% PURPOSE: Discrete element simulation of point particle
%  bouncing on a floor (at height 0) with
%  dissipation. The force between 
%  floor and particle is proportional to spring constant
%  k and the penetration into the floor (coordinate), 
%  as well as to the relative velocity to the floor.
%  The 
% USAGE: If the program is used as it is, the height is 
%  plotted at the time where the integrator choose to
%  evaluate the trajectory; if the lower part is commented
%  out, the solution will be interpolated ever at times which
%  are multiples of 0.1
% CAVEAT for the force computation: 
%  in  bouncing_ball_dissipation.m,
%  which is called by default, the possibility of non-smooth
%  variations in the damping force during approach 
%  are not avoided, so that for some solvers, the
%  maxima of the heights are unrealistic.
%  A function where the jump in the dissipative force
%  is prevented is bouncing_ball_dissipation7.4
% CAVEAT for the plotting:
%  With the current interpolation of the solution to
%  multiples of 0.05, not all points during the contact
%  with the floor will apper in the graphics;
%  If you use tspan=[0:10], all points will appear,
%  but while the particle is in free flight, only few
%  timesteps will be evaluated, so to obtain a 
%  continuous trajectory, it will be better to plot also
%  the trajectory with tspan=[0:.05:10] in the same graph
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

clear, format compact
global k, k=100;
global g, g=-9.81;
global D, D=0.1
x0=4
v0=1
tspan=[0:.05:10];
[t,y]=ode23('bouncing_ball_dissipation',tspan,[v0 x0]);
% Use this to limit the damping force during approach
%[t,y]=ode23('bouncing_ball_dissipation7_4',tspan,[v0 x0]);
plot(t,y(:,2),'k*-')

return
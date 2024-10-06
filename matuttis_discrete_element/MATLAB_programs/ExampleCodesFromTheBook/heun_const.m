function [t_plot,y_plot]=heun_const(fname,tspan,y0,dt);
% PURPOSE: Heun-Method (Runge-Kutta 2. Order)
%  with constant timestep
% USAGE: Call from ch3p148frictionaloscillator
%  The first three arguments are the same as for
%  MATLAB's built in integrators ode23, ode45 ....,
%  but the time-step must be specified as 4th argument
% CAVEAT: The timestep limits the resolution so in the 
%  worst case, the onset of sticking or slipping may be 
%  off by  nearly a full timestep. Nevertheless, as for 
%  these phenomena the experimental accuracy of the 
%  friction coefficient is the most crucial limitation 
%  of the accuracy, the limitation due to the timestep 
%  is not of practical consequence.
% LITERATURE: sec 2.3 and Tab. 2.1
% REVISION HISTORY: 26-May-2014 H.-G. Matuttis

global particle_sticks, 
       particle_sticks=0

if (tspan(1)>tspan(2))
  error('first entry of tspan must be smaller than second entry')
end

% Recompute the timestep to obtain an integer number
% of constant timesteps for the time interval which
% is close to the inputted timestep
dtold=dt;

nstep=round((tspan(2)-tspan(1))/dt);
dt=(tspan(2)-tspan(1))/nstep;
if (abs(dtold-dt)>1e-9*dt)
  warning(['dt=' num2str(dtold,8) ' has been changed to dt=' num2str(dt,8) ' to obtain the same timestep everywhere' ])  
end

y=y0;
ny=length(y);
t_plot=zeros(1,nstep);
y_plot=zeros(length(y),nstep);
t=tspan(1);
for istep=1:nstep
  t_plot(istep) = t;         % save time and y-values
  y_plot(1:ny,istep) = y;   % for plotting
  t_0=t;
  y_0=y;      
  k_0=feval(fname,t_0,y_0);
  dt_2=dt/2;
  y_1=y_0+dt_2*k_0;
  t_1=t_0+dt_2;
  k_1=feval(fname,t_1,y_1);
  y_2=y_0+dt*k_1;
  t_2=t_0+dt;
  k_2=feval(fname,t_2,y_2);
  
  t=t+dt;
% Enforce constraint for static friciton  
%   y=y_2;
  y=solout(y_0,k_0,y_2,k_2,dt);
  t=t_0+dt;
end

return

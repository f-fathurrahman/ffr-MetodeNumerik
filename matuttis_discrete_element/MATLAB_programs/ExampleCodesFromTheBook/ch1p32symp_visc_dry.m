% PURPOSE: Computes numerically the trajectories and 
%  velocities of the linear oscillator with no, viscous 
%  and dry friction. The graphs plotted are equivalent
%  to Fig. 1.15, 1.16 on p.32.
% USAGE: Calls the subroutine linosc_symp_visc_dry
%  If you don't know how to use numerical integrators
%  in MATLAB, study appAp419driverode.m first.
% CAVEATS: 
%  1. The dry friction case is treated "by
%   rejection" (not "numerically exact" as explained
%   in chapter 3), i.e. when the velocity should be zero,
%   the friction force oscillates between +mu and -mu.
%   The timestep is considerably reduced, the the
%   computation for dry friciton case takes considerable
%   time beyond t=22 with the current initial conditions.
%   If e.g. the initial velocity is reduced, the velocity
%   will become zero earlier and the time consumption
%   will increase. 
%  2. The solutions are numerical, not exact.
%   In particular, for zero friction and dissipation, the
%   maximal amplitude will not be constant, and
%   for the case of dry friction
% TODO: Try out different integrators (for ode45, ode23t,
%  ode23s and ode15 only the function names must be changed). 
%  Observe how the timesteps change for different integrators.
%  For some integrators, the simulation may not
%  finish and terminate with an error message!
% REVISION HISTORY: 23-May-2014 H.-G. Matuttis

clear;  
format compact
global mass, mass=1;  % mass set for all three cases 
global k, k=1; % spring constant set for all three cases
global D   % damping, different for each case
global mu  % friction coefficient, different for each case

tspan=[0 22]
ystart=[0 2]';

D=0.;
mu= 0.0;
[t_symp,y_symp]=ode23('linosc_symp_visc_dry',tspan,ystart);
x_symp=y_symp(:,2);
v_symp=y_symp(:,1);

D      = 0.1;
mu     = 0.0;
[t_visc,y_visc]=ode23('linosc_symp_visc_dry',tspan,ystart);
x_visc=y_visc(:,2);
v_visc=y_visc(:,1);

warning('If mu is too large (or the initial velocity too small), the following computation may take quite long!')

D      = 0.0;
mu     = 0.15;
[t_dry,y_dry]=ode23('linosc_symp_visc_dry',tspan,ystart);
x_dry=y_dry(:,2);
v_dry=y_dry(:,1);

clf
subplot(2,1,1), hold on
plot(t_symp,x_symp,'k+--')
plot(t_visc,x_visc,'k*-','Color',[0.5 0.5 0.5])
plot(t_dry,x_dry,'kd-')
plot([0 21.4],[2 0],'k:')
plot([0 21.4],-[2 0],'k:')
plot([0 25],-[0 0],'k:')
     axis([0 25 -2 2])
     xlabel('time t')
     ylabel('x(t)')
legend('\mu=0,\delta=0',...
       '\mu=0,\delta=0.1',...
       '\mu=0.15,\delta=0')

subplot(2,1,2), hold on
plot(t_symp,v_symp,'k+--')
plot(t_visc,v_visc,'k*-','Color',[0.5 0.5 0.5])
plot(t_dry,v_dry,'kd-')
plot([0 21.4],[2 0],'k:')
plot([0 21.4],-[2 0],'k:')
plot([0 25],-[0 0],'k:')
axis([0 25 -2 2])
xlabel('time t')
ylabel(' v(t)')
legend('\mu=0,\delta=0',...
       '\mu=0,\delta=0.1',...
       '\mu=0.15,\delta=0')
  
return


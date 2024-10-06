% PURPOSE: Compute the trajectory and velocity of the 
%  stick-slip oscillator secion 3.3.2 "numerically 
%  exact" by detection of the static friction regime
%  via the flow of the phase space
% USAGE: calls stick_slip_oscillator.m,
%  integrates the ODE via heun_const.,
%  and implements static friction via
%  solout.m 
% CAVEAT: This problem uses an integrator
%  with constant timestep (Runge-Kutta Heun)
%  so that the timestep must be specified:
%  for too large timesteps, the solution
%  becomes meaningless
% LITERATURE: sec. 3.3.2
% REVISION HISTORY: 26-May-2014 H.-G. Matuttis

clear all;  
format compact
global particle_sticks

global mass,   mass=1;
global spring, spring=1;
global D,      D     =0.1;
global my,     my    =4; 
global A,      A     =1;   
global omega,  omega =pi;

disp('    ')
disp('"Numerically exact" solution of a stick-slip oscillator problem')
dt=.05 % Timestep

ystart=[4 3]'; % Initial condition: v(0)=4,x(0)=3
hold on     
[t_plot,y_plot]=heun_const('stick_slip_oscillator',[0 10],ystart,dt);
clf 
% Plot velocities
plot(t_plot,y_plot(1,:),'+-','Color',[.2 .2 .2]);
hold on
% Plot trajectories into the same plot
plot(t_plot,y_plot(2,:),'x-','Color',0*[.2 .2 .2]);
legend(' x(t) ',' v(t) ')
xlabel('time')

return
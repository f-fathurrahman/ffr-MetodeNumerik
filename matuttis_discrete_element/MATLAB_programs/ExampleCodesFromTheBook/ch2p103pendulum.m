% PURPOSE: Compute and plot the trajectories
%  for a pendulum implemented via constraints
%  (differential algebraic equations, DAE)
%  with consistent and inconsistent initial
%  conditions
% USAGE: Calls the function pendulum.m which
%  contains the actual equation.
% CAVEAT: If the initial conditions for the
%  velocity deviate too strongly from 
%  "consistent" initial conditions, the 
%  solution will spiral to infinity very
%  fast, and the graph will look meaningless
% TODO: For consistent initial conditions,
%  replace the usage of ode23 with ode45 and
%  observe hove the trajectory is maintained
%  with higher accuracy. Further, modify
%  the relative and absolute tolerance
%  (see appAp419driverode.m for an example)
%  and see how the accuracy for the constraint
%  improves.
% REVISION HISTORY: 23-May-2014 H.-G. Matuttis

clear
format compact

disp('Computation of a pendulum trajectory in DAE-formulation')
global m, m=10; % mass
global g, g=10; % gravitation
tspan=[0 10]; % Time interval for all simulations

% consistent initial condition for position
% (length l=1)
% and velocity, trajectory of a circular arc
intitialcondition=[1 0 0 0]; 
[t,y] = ode23('pendulum',tspan,intitialcondition,odeset('reltol',1e-5));
clf
subplot(2,1,1)
plot(y(:,1),y(:,2),'--','Color',[0 0 0])
text(1,-.5,'l=1')

% consistent initial condition for velocity
% and inconsistent initial condition for
% position (length l=2), which nevertheles has 
% no effect, as the position is not used in
% the implementation of the constraints:
% trajectory of a circular arc
intitialcondition=[2 0 0 0];
[t,y] = ode23('pendulum',tspan,intitialcondition,odeset('reltol',1e-5));
hold on
plot(y(:,1),y(:,2),'-','Color',.3+[0 0 0])
text(1.01,-1.8,'l=2')
axis equal

% consistent initial condition for position
% (length l=1) and inconsistent initial 
% conditions for the velocity, trajectory
% spirals away from the correct one for 
% a circular arc
intitialcondition=[1 0  .1  1];
[t,y] = ode23('pendulum',tspan,intitialcondition,odeset('reltol',1e-5));
subplot(2,1,2) 
plot(y(:,1),y(:,2),'-','Color',[0 0 0])
text(1.01,-1.8,'l=2')

% consistent initial condition for position
% (length l=1) and inconsistent initial 
% conditions for the velocity with
% larger error than before, trajectory
% spirals away from the correct one for 
% a circular arc
intitialcondition=[1 0  1  1];
[t,y] = ode23('pendulum',tspan,intitialcondition,odeset('reltol',1e-5));
hold on
plot(y(:,1),y(:,2),'--','Color',.3+[0 0 0])
text(1.01,-1.8,'l=2')

axis equal



return




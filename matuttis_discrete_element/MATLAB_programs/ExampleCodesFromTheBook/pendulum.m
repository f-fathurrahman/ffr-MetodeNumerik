function [dydt]=pendulum(t,y);
% PURPOSE: Ordinary differential equation
%  for the trajectory of a pendulum implemented 
%  with constraints
% USAGE: Call from ch3p103pendulum.m
% IMPLEMENTATION: The coordinates and 
%  velocities in the y-vector are entered as
%   y(1)=r_x;
%   y(2)=r_y;
%   y(3)=v_x;
%   y(4)=v_y;
% CAVEAT: The stabilization by projection can
%  not be implemented in this file, and it 
%  cannot be implemented with MATLAB's 
%  ode-solvers; One needs (has to write on one's own)
%  a solver where one has access to the  function
%  value computed at the end of a timestep. 
% REVISION HISTORY: 23-May-2014 H.-G. Matuttis
global m
global g

%f=[0 -m*g];

% Computation of the Lagrange multiplier
lambda=(-f*[y(1) y(2)]'-m*[y(3) y(4)]*[y(3) y(4)]')/([y(1) y(2)]*[y(1) y(2)]');

dydt(1,1)=y(3);
dydt(2,1)=y(4);
dydt(3,1)=lambda*y(1)/m;
dydt(4,1)=(lambda*y(2)-m*g)/m;

return
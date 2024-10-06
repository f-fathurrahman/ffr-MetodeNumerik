function [dy]=bouncing_ball(t,y)
% PURPOSE: function for bouncing ball without dissipation,
%  the force acts if the particle is below 0,
%  the force is proportional to the distance to 0 and
%  to the spring constant k.
% USAGE: Use with ch7p225roundparticledem7_1
% Program 7.2 in the book
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

global g
global k
if (y(2)>=0)
  dy=[g
  y(1)];
else
  dy=[g-k*y(2)
  y(1)];
end
return

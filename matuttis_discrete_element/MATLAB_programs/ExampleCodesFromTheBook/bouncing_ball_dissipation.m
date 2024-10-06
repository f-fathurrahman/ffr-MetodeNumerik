function [dy]=bouncing_ball_dissipation(t,y)
% PURPOSE: function bouncing ball with dissipation
% ALGORITHM: the force acts if the particle is below 0,
%  the elastic force is proportional to the distance to 0 
%  and to the spring constant k,
%  the damping force is proportional to D and the velocity.
%  Additionally, a condition is used to cut of the
%  damping force if it overcompensates the elastic force.
% USAGE: Use with ch7p225roundparticledem7_3
% Program 7.4 in the book
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

global g
global k
global D



if (y(2)>=0)
  dy=[g
  y(1)];
else
  f_el=-k*y(2);
  f_damp=-D*y(1);
  f_tot=f_el+f_damp;
% Comment this out to obtain unphysical solutions
% where due to large damping force, the interaction
% can become atractive
  if (sign(f_tot*f_el)<0)
    f_tot=0;
  end
  dy=[g+f_tot
  y(1)];
end
return

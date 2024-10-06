function Y=solout(Y0,F0,Y1,F1,dt);
% PURPOSE: Project the solution onto the constraint manifold 
%  for static friction, i.e. set the velocity to zero and
%  the position to the position at the beginning of the 
%  timestep if the conditions for static friction are met
%  see discussion sec. 3.3.1
% USAGE: To be called from an integrator at the end
%  of the timestep
%  Y0 is the solution at the beginning of the timestep,
%  F0 its Jacobean (gradients),
%  Y1 is the solution at the end of the timestep,
%  F1 its Jacobean (gradients)
%  and dt the timestep
%  The components are
%    Y..(1) velocity, F..(1) acceleration
%    Y..(2) position, F..(2) velocity
% ALGORITHM: 
%  0. Set Y as the Heun-solution.
%     Replace it if sticking is detected.
%  1. If sticking has been determined in the previous
%   timestep, set a flag for the next timestep
%  2. If sticking has been determined for the current
%   timestep, set a flag for the next timestep
%   If the sticking flag is set, set the position to the 
%   position at the beginning of the current timestep
%   and the velocity to 0
% CAVEAT:
%  While the flow condition for the Painleve Paradox
%  (aI>0)&(aII>0) will not occur for interpolation
%  between known values, it may occur for extrapolation
%  with low accuracy (i.e. the numerical error and the
%  inaccuracy of the extrapolation, due to the non-smooth
%  evolution of the force is so large that the sign
%  in both aI or aII becomes unreliable).
%  Therefore, the second condition uses
%     if (aI*aII<0)&~(abs(vnew)<1e-4)
%  instead of 
%     if ((aI<0)&(aII<0)) 
% Extrapolate the solution for the next time
%  step, Determine whether there is a velocity reversal
%  solution goes through 
%  0 for the velocity, and if 
% Extrapolation: Predict wether the next step will go through 0
% REVISION HISTORY: 26-May-2014 H.-G. Matuttis


global particle_sticks

YOLD=Y0;
Y=Y1; 
if (particle_sticks==1)
  Y(1)=0;  
  Y(2)=YOLD(2);
end
delta=1e-4;
aI = (Y0(1)+delta*F0(1))/delta;
aII=-(Y1(1)+delta*F1(1))/delta;
% Condition for static friction
if ((aI<0)&(aII<0)) 
  Y(1)=0;
  Y(2)=YOLD(2);
  particle_sticks=1;
else
  particle_sticks=0;
end
    
% Predict the velocity in the next time-step using
% the dense output (Taylor-series using the initial
% and final values y0, y1, as well as the initial and 
% final gradients f0 and f1)
y0=Y0(1);
f0=F0(1);
y1=Y1(1);  
f1=F1(1);
y1my0=y1-y0;
% Extrapolated Velocity for the next timestep
vnew=y0+(dt*f0+(3*y1my0-dt*(2*f0+f1)+(-2*y1my0+dt*(f0+f1))*2)*2)*2;
if (vnew*y1<=0) % Velocity reversal in the next interval
  delta=1e-4;
% Extrapolated acceleration for the ntext timestep
  anew=F1(1)+dt*(F1(1)-F0(1));
  aI = (Y1(1)+delta*F1(1))/delta;
  aII=-( vnew+delta*anew)/delta;
% Set sliding friction if the conditions for
% sliding friction aI*aII<0 are met, but
% only if the predicted velocity is not
% numerically zero (in this case, larger than 10e-4
  if (aI*aII<0)&~(abs(vnew)<1e-4)
    particle_sticks=0;
  else
    Y(1)=0;
    Y(2)=YOLD(2);
    particle_sticks=1;
  end
else
  particle_sticks=0;
end

return
    
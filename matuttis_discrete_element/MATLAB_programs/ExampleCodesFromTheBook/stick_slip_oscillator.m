% Model:
% Particles are linked with a spring of linear force
% Numbering starts at the left !
% Convex hull version, see sec. 3.3.1
% ALGORITHM: First compute the Velocity-Independent forces
% REVISION HISTORY: 26-May-2014 H.-G. Matuttis

function [dy] = stick_slip_oscillator(t,y)
% Return dy/dt = f(t,y).
global mass   
global spring 
global D     
global my     
global A      
global omega  



v=y(1);
% the solution is quite robust with
% respect to modifications of this delta:

force=A*cos(omega*t)-2*D*y(1)-spring*y(2);
forceI = force - my;
forceII =force + my;     
delta=1e-4;
aI = (v+delta*forceI/mass)/delta;
aII=-(v+delta*(forceII/mass))/delta; 
if (aI*aII<0) % Flow traverses: dynamic friction
  if (v>0)
    force=forceI;
  else
    force=forceII;
  end
elseif ((aI>0)&(aII>0)) % Painleve Paradox
  error('Flow pulls to both sides')
elseif ((aI<0)&(aII<0)) % Static friction
  lambda=aI/(aI+aII);
% Solution "in the manifold"  
  force=(1-lambda)*forceI+lambda*forceII;
  v=0; 
end
dy(1,1) = force/mass;
dy(2,1) = v ;

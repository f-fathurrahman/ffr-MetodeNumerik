% PURPOSE: Plot the solution of the ordinary 
%  differential equation y'=y evaluated
%  at the default points choosen by the integrator,
%  and interpolated to multiples of 0.1 of the
%  abscissa values
% USAGE: Driver for the function expfun, which
%  numerically integrates out y'=y
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

clear
format short
format compact
xinitial=1
% time steps between 0 and 3 will be selected 
% by the integrator:
timespan1=[0 3] 
% time steps between 0 and 3 will be selected 
% by the integrator, but the solution will 
% be interpolated to the multiples of 0.1,
% and only these values will be outputted:
timespan2=[0:0.1:3];
[t1,x1]=ode23('expfun',timespan1,xinitial);
[t2,x2]=ode23('expfun',timespan2,xinitial);
clf
plot(t1,x1,'*',t2,x2,'o')
legend('computed solution','interpolated solution')

return
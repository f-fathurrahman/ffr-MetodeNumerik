% PURPOSE: Plot the solution of the ordinary 
%  differential equation y'=y with different
%  accuracies
% USAGE: Driver for the function expfun, which
%  numerically integrates out y'=y
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

clear
format short
format compact
xinitial=1, timespan=[0 3]
options1=odeset('RelTol',1e-1,'AbsTol',1e-1);
[t1,x1]=ode23('expfun',timespan,xinitial,options1);
options2=odeset('RelTol',1e-4,'AbsTol',1e-4);
[t2,x2]=ode23('expfun',timespan,xinitial,options2);

clf
plot(t1,x1,'*--',t2,x2,'o-')
legend('RelTol=10^{-1},AbsTol=10^{-1}',...
       'RelTol=10^{-4},AbsTol=10^{-4}')


return
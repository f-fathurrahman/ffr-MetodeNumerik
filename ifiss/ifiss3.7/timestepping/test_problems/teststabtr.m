%TESTSTABTR sample STABTR run
%   X-IFISS scriptfile: DJS; 15 March 2018
% Copyright (c) 2018 D.J. Silvester, Huining Yang
fprintf('Generating steady state solution... ')
batchmode('CD4')
load batchrun
fprintf('\n\nTimestepping to quiescent state ... '), pause
close all
nb=length(f);
[DT,U,Udot,time] = stabtrX(Asupg,Q,f,xsupg,1e-9,100,5e-5,10,1);
figure(2), loglog(time,DT,'bx'), axis square
xlabel('time'), ylabel('timestep'), title('time step evolution')
pause
tsolmovie(U,time,xy,x,y,19,'jet');


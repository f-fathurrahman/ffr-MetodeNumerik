%RESTART_STEP_NAVIER restart unsteady flow in step domain
%   s-IFISS scriptfile: DJS; 2 August 2023.
% Copyright (c) 2014 D.J. Silvester, H.C. Elman
% Modified (c) 2023 M.D. Mihajlovic; corrected DJS 9 November 2023
clear all; close all
global pde viscosity

fprintf('Navier-Stokes flow in step domain ... \n')
fprintf('Restarting integration from Checkpoint datafile ...\n')
% load data from original integration
gohome, cd datafiles
load step_unsteadyflow.mat
load stabtrNS_end.mat
fprintf('viscosity parameter is  %11.3e\n',viscosity)
restart=initdata.restart+1;
dt=initdata.dt; tbegin=initdata.time;
n=length(dt)-1; oldDT=DT;
dt0=initdata.dt0;
u=initdata.u;

% load assembled matrices
gohome, cd datafiles
load step_stokes_nobc.mat
% initialize for time stepping iteration
fprintf('restarting from %11.3e seconds\n',tbegin)
%% set parameters
pde=14; domain=1;
tfinal = default('new final time? (1e14)',1e14);
nmax = default('number of timesteps (x200)? (default 5)',5); %% 5 | 1000 timesteps
tol = default('accuracy tolerance? (default 3e-5)',3e-5);
nonlin = default('number of Picard steps? (default 1)',1);
nstar = default('averaging frequency? (default 10)',10);
vswitch = default('plot solution evolution? 1/0',0);
if(vswitch==1)
   vfreq=default('plotting frequency (default every 10 time steps)',10);
   fps=default('movie frame rate [fps] (default 16)',16);
   xc=default('x coordinate for the horizontal velocity profile (default L/2)',outbnd/2);
   if(xc<0 | xc>outbnd), error('Oops.. x coordinate outside of the domain'), end
   xr=x-xc; [xmin,ix]=min(abs(xr)); xref=x(ix);
end
gzero=zeros(size(g));
%% compute solution
tstart = tic;
np=length(g); 
AxB='defaultAxB';
[DT,U,Udot,xtime] = stabtrNS(qmethod,xy,mv,bound,A,B,sparse(np,np),G,AxB,...
                     u,dt0,tfinal,tol,nstar,1,nonlin,nmax,restart);
etoc=toc(tstart); fprintf('Integration took  %8.3e seconds\n\n',etoc) 
save step_unsteadyflow.mat DT U Udot time viscosity
clear U Udot
% visualise the timestep sequence plot
marker = 'k.';
dttplot(oldDT(1:end-1),99,marker,0)
marker = 'r.'; 
dttplot(DT,99,marker,length(oldDT)-1)
pause(2)
%
% postprocessing
if(vswitch==1)
stepflowmovie
end
return

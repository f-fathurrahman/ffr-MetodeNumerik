%UNSTEADY_STEP_NAVIER solve Navier-Stokes problem in step domain
%   IFISS scriptfile: DJS; 22 August 2016, AR; 09 January 2019
% Copyright (c) 2009 D.J. Silvester
% Modified (c) 2023 M.D. Mihajlovic; corrected DJS 26 December 2024
%clear variables
global pde
fprintf('Unsteady flow in a step domain ...\n')
global viscosity
viscosity=default('viscosity parameter (default 1/220)',1/220);
global DELTA
DELTA=default('inlet perturbation magnitude (default 0)',0);
%
%% initialize for time stepping iteration: compute Stokes solution
%% load assembled matrices
gohome
cd datafiles
load step_stokes_nobc.mat
%
%
%% set parameters
pde=14; domain=3;
fprintf('Discrete Saddle-Point DAE system ...\n')
tfinal = default('target time? (default 1e8)',1e8);
nmax = default('number of timesteps (x200)? (default 5)',5); %% 5 | 1000 timesteps
tol = default('accuracy tolerance? (default 3e-5)',3e-5);
nonlin = default('number of Picard steps? (default 2)',2);
nstar = default('averaging frequency? (default 10)',10);
vswitch = default('plot flow solution evolution? 1/0',0);
if(vswitch==1)
   vfreq=default('plotting frequency (default every 10 time steps)',10);
   fps=default('movie frame rate [fps] (default 16)',16);
   xc=default('x coordinate for the horizontal velocity profile (default L/2)',outbnd/2);
   if(xc<0 | xc>outbnd), error('Oops.. x coordinate outside of the domain'), end
   xr=x-xc; [xmin,ix]=min(abs(xr)); xref=x(ix);
end
dtzero=1e-9; uzero=zeros(size(f));gzero=zeros(size(g));
%% compute solution
tstart = tic;
if qmethod>1,
np=length(g);
AxB='defaultAxB';
[DT,U,Udot,time] = stabtrNS(qmethod,xy,mv,bound,A,B,sparse(np,np),G,AxB,...
                            uzero,dtzero,tfinal,tol,nstar,1,nonlin,nmax,0);
else
error('Oops.. functionality is not available for stabilized approximations')
end
etoc=toc(tstart); fprintf('Integration took  %8.3e seconds\n\n',etoc)
save step_unsteadyflow.mat DT U Udot time viscosity
clear U Udot
%
% visualise the timestep sequence plot
marker = 'k.';
offset=0;
dttplot(DT(1:end-1),97,marker,offset)
pause(1)
%
%  postprocessing
if(vswitch==1)
stepflowmovie
end
return

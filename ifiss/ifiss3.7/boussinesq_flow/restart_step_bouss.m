%RESTART_STEP_BOUSS re-solve Boussinesq flow problem in step (PC version)
%   IFISS scriptfile: DJS; 13 May 2012.
% Copyright (c) 2012 D.J. Silvester, M.D. Mihajlovic.
% Modified (c) 2023 M.D. Mihajlovic
% Modified (c) 2023 M.D. Mihajlovic; corrected DJS 10 November 2023

clear variables
global Ra Pr H L 
%
% initialize for time stepping iteration: load assembled matrices
gohome; 
cd datafiles; 
load step_bouss_nobc.mat; load step_grid1h.mat
% load the restart data and check domain parameter
load stabtrBouss_restart.mat
if domain~=3;
error('Oops: the domain seems to be incompatible!'), end

fprintf('Unsteady Boussinesq flow in a backward step domain ...\n')
fprintf('Restarting integration from Checkpoint datafile ...\n')
Ra=default('Rayleigh number (default is unchanged)',Ra);
Pr=default('Prandtl number (default unchanged)',Pr);
fprintf('Ra is %g and Pr is %g \n\n',Ra,Pr); 

% save timestep and time history
oldtime=time; oldDT=DT; oldtol=tol;
clear U Udot T Tdot DT time solDT soltime

% set new parameters
fprintf('integration will be restarted from  %11.3e seconds\n',t)
tfinal = default('new target time? (default 200)',200);
nmax = default('number of timesteps (x200)? (default 5)',5); %% 5 | 1000 timesteps
tol = default('accuracy tolerance? (default 3e-5)',3e-5);
nonlin = default('number of nonlinear Picard correction steps? (default is none)',0);
nstar = default('averaging frequency? (default 10)',10);
vswitch = default('plot solution evolution? 1/0 (default is yes)',1);
if(vswitch==1)
   vfreq=default('plotting frequency (default every 10 time steps)',10);
   fps=default('movie frame rate [fps] (default 16)',16);
   xc=default('x coordinate for the horizontal velocity profile (default L/2)',L/2);
   if(xc<-1 | xc>L), error('Oops.. x coordinate outside of the domain!'), end
   xr=x-xc;                  %  find the closest grid line to xc
   [xmin,ix]=min(abs(xr));
   xref=x(ix);               %  find the points on the referent grid line
end
xout = default('generate solution history data points? 1/0 (default is no)',0);
tout = xout+1; 
%
% unpack grid data
xyv=grid(1).xyv;  [nnv,dd]=size(xyv(:,1));
xyp=grid(1).xyp;  [nnp,dd]=size(xyp(:,1));
xyt=grid(1).xyt;  [nnt,dd]=size(xyt(:,1));
uzero=zeros(2*nnv,1); gzero=zeros(nnp,1); hzero=zeros(nnt,1);
% 
% compute solution
tic,
if qmethod==12,
[refDT,U,Udot,T,Tdot,reftime] = ...
	stabtrBouss(nonlin,qmethod,grid,spmat,initdata,tfinal,tol,nstar,nmax,tout);
else
error('Oops.. mixed approximation strategy is not yet implemented!')
end
etoc=toc; fprintf('integration took  %8.3e seconds\n\n',etoc) 
%
% update the time history in the checkpoint file if it exists
if exist('stabtrBouss_restart.mat','file')==2,
xtime=load('stabtrBouss_restart.mat','time','DT');
time=[oldtime,xtime.time]; DT=[oldDT,xtime.DT];
save('stabtrBouss_restart.mat','time','DT','pde','domain','-append'); end
% update the problem definition in the final solution file
save('stabtrBouss_end.mat','time','DT','pde','domain','-append');
%
% update the timestep sequence and plot
fprintf('Timestep history: \noriginal ... ') 
DT=[oldDT,refDT]; time=[oldtime,reftime];
marker = 'kx'; offset=0;
dttplot(oldDT,199,marker,offset)
fprintf('final ... ') 
marker = 'rx';
dttplot(refDT,199,marker,length(oldDT)-1)
pause(2)
%
% postprocessing
if vswitch==1
stepboussmovie
end
return

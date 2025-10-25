function [DT,U,Udot,time] = stabtrNS(qmethod,xy,ev,bound,A,B,C,G,AxB,...
                                uzero,dtzero,tfinal,tol,nstar,tout,nonlin,nmax,restart)
%STABTRNS Navier-Stokes integrator using stabilized TR 
%   [DT,U,Udot,time] = stabtrNS(qmethod,xy,ev,bound,A,B,C,G,AxB,...
%                          uzero,dtzero,tfinal,tol,nstar,0,nonlin,nmax,restart);
%   input
%          qmethod   approximation method
%          xy        vertex coordinate vector
%          ev        mv/ev  Q2/Q1 element mapping matrix
%          bound     boundary vertex vector 
%          A, B, C   matrices defining the saddle point system
%          G         mass matrix for velocity
%          AxB       saddle point system solver function_name
%          uzero     initial velocity condition
%          dtzero    initial timestep
%          tfinal    final time 
%          tol       local accuracy error tolerance   
%          nstar     averaging frequency   
%          tout      output data switch
%          nonlin    number of nonlinear corrections (optional, default is 0)
%          nmax      number of integration intervals (*nnt)
%          restart   restart flag
%   output
%          DT        timestep history
%          U         solution history
%          Udot      solution time derivative history
%          time      discrete time evolution
%
% calls functions: unsteadyflowbc, dtflowbc, resetflowbc, AxBhandle
%   IFISS function: DJS; 7 June 2010; ... 17 August 2016 .. 27 December 2024
% Copyright (c) 2006 D.J. Silvester, D.F. Griffiths, M. Schneider
% Modified (c) M.D. Mihajlovic; corrected DJS 9 November 2023
%
if nargin < 16, nonlin=0; end
global viscosity
%%% set the restart parameter
nnt=200; 
%%%
warning('on','MATLAB:nearlySingularMatrix'); wchar=27;
fprintf('Solving discrete N-S system using stabilized TR ...\n')
AxBhandle=str2func(AxB),
[np,nuv]=size(B); nv=nuv/2;
ub=uzero; 
Ndof=length(ub); T=tfinal; dt0=dtzero;
% preallocate array dimensions
nntmax=nmax*nnt;
DT = zeros(1,nntmax); time = zeros(1,nntmax); 
U = zeros(Ndof,nnt); Udot = zeros(Ndof,nnt);
%
% initialise file used for movie generation
if(exist('VelocityScaling.mat')==2), delete VelocityScaling.mat, end
uxmax=-inf; uymax=-inf; uxmin=inf; uymin=inf;
gohome, cd datafiles
save('VelocityScaling.mat','uxmax','uymax','uxmin','uymin');


% pack start data
krestart=0;
uxmax=0; uymax=0;
%%% zeroth time step 
f=zeros(nuv,1);gzero=zeros(np,1);    
%% unpack initialization structure
if(restart==0)
   fprintf('intermediate (CheckPoint) datafile: stabtrNS_end\n')
   initdata=struct('u',uzero,'ub',ub,'udot',Udot,'udotb',Udot, ...
                   'dt0',dtzero,'dt',dtzero,'time',time,'n',0,'restart',restart); 
   ub=initdata.ub; dt0=initdata.dt0; % dt0=dtzero;
   n=0;
%
% initialization: potential flow solve
   [fst,gst] = unsteadyflowbc(viscosity*A+(1/dt0)*G,B,f,gzero,xy,bound,dtzero);
   bdy=fst(bound);
   if(tout==1)
      fprintf('\n  initial nonlinear residual is %e \n',norm([fst;gst]));
      fprintf('             boundary change is %e \n',norm(bdy));
   end
   if(qmethod>1)  
      N = navier_q2(xy,ev,ub,1);
   elseif(qmethod<=1) 
      N = navier_q1(xy,ev,ub,1); 
   end
   Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Gst,Bst,fzz,gzz] = dtflowbc(G,B,-Anst*ub,gzero,xy,bound,0,dt0);
   xpot = [Gst,Bst';Bst,-C]\[fzz;gzz];
   udotb=xpot(1:nuv); pdotb=xpot(nuv+1:end); 
% compatibility check
   divres_norm=norm(gzero-Bst*ub,inf);
   if divres_norm>10*eps
      fprintf('initial divergence residual is %e',divres_norm)
      error('incompatible  initial data!')
   end
% take care of system singularity warnings
   [lastmsg, msgid] = lastwarn; 
   if(length(msgid)==wchar)
      lastmsg, fprintf('This should not cause difficulty for enclosed flow problems.\n');
   end
   warning('off','MATLAB:nearlySingularMatrix'); 
%
%%% first time step
   n=1; k=1;                                         % time step index 
   ww = ub+dt0*udotb;                                % w(dt0)
   if(qmethod>1)  
      N = navier_q2(xy,ev,ww,0);
   elseif(qmethod<=1) 
      N = navier_q1(xy,ev,ww,0); 
   end
   Anst = 2*G + dt0*viscosity*A + dt0*[N, sparse(nv,nv); sparse(nv,nv), N];
   fnst = G*udotb - (viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N])*ub;
   gnst= -(B*ub); %fprintf('\ndivergence residual is %e',norm(gnst))
   [Anst,Bst,fzz,gzz] = dtflowbc(Anst,B,fnst,(1/dt0)*gnst,xy,bound,0,dt0);
%---------- compute unscaled pressure using AxB
   [v,pns] = AxBhandle(Anst,Bst,C,fzz,gzz,1);
   xns=[v;pns];                                     % first TR step
   u = ub + dt0*v;                                  % u(dt0) 
   udot = 2*v - udotb;	                             % du/dt(dt0)
   udd = (udot-udotb)/dt0;		                     % second derivative
   dt = dt0;    t = dt;
   acc = sqrt((udot'*(G*udot)));
   normrnst=0; % --- set initial momentum residual to zero
   inner=setdiff([1:nuv]',bound); % --- interior velocity nodes
   if(nonlin==0)
      fprintf('\n StabTR with no nonlinear corrections');
   else
      fprintf('\n StabTR with %2i nonlinear corrections',nonlin); 
   end
   if(tout==1) 
      fprintf('\n   step  timestep       time        residual      divresidual     acceleration'),
      fprintf('\n  %4i   %9.3e   %11.3e    %9.3e     %9.3e      %11.3e', n, dt, t+dt,normrnst, norm(gnst),acc);
   end
%
%%% second time step
   n=2; k=2;
   DT(1:2) = [dt0, dt0];	
   U(:,1) = ww; U(:,2) = u;
   Udot(:,1) = udotb; Udot(:,2) = udot;
   time(1:2) = [0,dt0];
%------------   n=n+1; k=k+1;
else   %%%% hot restart
   load stabtrNS_end.mat
   n=initdata.n; t=initdata.time(end);
   tbegin=initdata.time; dt0=initdata.dt0; dt=initdata.dt;
   u=initdata.u; ub=initdata.ub;
   udot=initdata.udot; udotb=initdata.udotb;
   fprintf('\n\nTimestep %5i  Time %11.3e \n',n,t)
   n=1; k=1;
   udd = (udot-udotb)/dt0; %  tdd = (tdot-tdotb)/dt0;	
   time(1) = t;  DT(1) = dt;	
   U(:,1)  = u;  Udot(:,1) = udot; 
   n=n+1; k=k+1;
   fprintf('Restarting integration ...\n ')
   if(nonlin==0)
      fprintf('\n StabTR with no nonlinear corrections');
   else
      fprintf('\n StabTR with %2i nonlinear corrections',nonlin); 
   end
   if(tout==1)
    fprintf('\n   step  timestep       time        residual      divresidual     acceleration'),
                 end
end

%---------------- loop until time limit is reached
flag = 0; nrej=0; avflag = 0;  nav=nstar; tstar=1;
avstep=0;
while t <= T  & flag==0   
   if(t+dt>T) 
      dt = T-t; 
      flag = 1; 
   end          % fix final time step
   if(n==nntmax) 
      flag=1; 
      fprintf('\nToo slow -- step limit reached!')
   end
%%%%%% general TR step
   ww = (1+(dt/dt0))*u - (dt/dt0)*ub;
   if(qmethod>1)  
      N = navier_q2(xy,ev,ww,0);
   elseif(qmethod<=1) 
      N = navier_q1(xy,ev,ww,0); 
   end
   Anst = 2*G + dt*viscosity*A + dt*[N, sparse(nv,nv); sparse(nv,nv), N];
   fnst = G*udot -(viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N])*u;
   gnst= -(B*u);
   %  if t>0,gnst=gzero;end
%-------- correction code
   dtz=dt;
   if avstep==1, dtz=dt-dtx; end
   [Anst,Bst,fzz,gzz] = dtflowbc(Anst,B,fnst,(1/dtz)*gnst,xy,bound,t,t+dtz);
%---------- compute scaled pressure using AxB
   [v,pns] = AxBhandle(Anst,Bst,C,fzz,gzz,n); 
   xns=[v;pns];                               % general TR step
%
%%%%% nonlinear iterations
   if(nonlin>0)
      for itpic=1:nonlin    % perform Picard iteration(s)
         iuTR  = u + dt*v;
         N = navier_q2(xy,ev,iuTR,0);
         Anst = 2*G + dt*viscosity*A + dt*[N, sparse(nv,nv); sparse(nv,nv), N];
         fnst = G*udot -(viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N])*u;
         gnst= -(B*u);
         [Anst,Bst,fzz,gzz] = dtflowbc(Anst,B,fnst,(1/dtz)*gnst,xy,bound,t,t+dtz);
         %---------- compute scaled pressure using AxB
         [vpic,pns] = AxBhandle(Anst,Bst,C,fzz,gzz,n);
         fprintf('\n  %8.2e   --- nonlinear correction  ', norm(vpic-v));
         v=vpic;
      end
   end
   w = udot + (.5*dt)*udd;
   udiff = v - w;
   uAB2  = u + dt*w;
   upTR  = u + dt*v;
% local truncation error estimate   
   d = (dt^2/(3*(dt+dt0)))*sqrt((udiff'*(G*udiff)));
% acceleration
   acc = sqrt((udot'*(G*udot)));
%%%%%% time step rejection
   if(d < ((1/0.7)^3)*tol)
%     accepted step
       if(t>tstar & ~rem(n,nav))    % smooth by averaging
         ub  = .5*(u+ub);
%---     ub  = resetflowbc(ub,xy,bound,t-.5*dt0); % don't reset boundary values
         udotb = .5*(udot+udotb);
	     dt0 = .5*(dt+dt0);                     %% correct the old timestep
         u = .5*(u + upTR);
         udot = v;
         t = t + .5*dt;      dtx=0.5*dt;             % leave dt unchanged
%--- correction  u  = resetflowbc(u,xy,bound,t);     % reset boundary value
         avflag=1; avstep=1;
         if(nav == 1e4) 
            nav = n; 
         end
         if(tout==1)
            fprintf('\n  %4i   %9.3e   %11.3e  ', n, .5*dt, t);
            fprintf('--- Averaging step');
         end
      else
%     regular step
         avstep=0;
         dt0 = dt;
         t = t+dt0;
         ub = u;
         u = upTR;
         udotb = udot;
         udot = 2*v - udot;
      end
      udd = (udot-udotb)/dt0;
% --- residual calculation
      rnst = G*udot + (viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N])*u + B'*pns;
      rnst(bound)=0; rnst(nv+bound)=0;
      gnst= -(B*u);
      if(tout==1)
         fprintf('\n  %4i   %9.3e   %11.3e    %9.3e     %9.3e      %11.3e', ...
                      n, dt, t,norm(rnst), norm(gnst),acc);
      end

%%%%%% save solution data
      DT(n) = dt; U(:,k) = u;  time(n) = t;
      Udot(:,k) = udot;
      if(avstep)    %    correct old values after averaging step
         U(:,k-1) = ub; Udot(:,k-1) = udotb; 
      end
      if(mod(k,nnt)==0)       %  need to reorder memory
         timeptr=restart*nnt+1;
         krestart=krestart+1; dateflag=date; restart=restart+1;
         initdata=struct('u',u,'ub',ub,'udot',udot,'udotb',udotb,'dt0',dt,...
                         'dt',dt*(tol/d)^(1/3),'time',t,'n',n-1,'restart',restart');
         gohome; cd datafiles
         longtime=time; longDT=DT;
         kknz=[1,find(time)]; time=time(kknz); DT=DT(kknz);
      %  soltime=time(end-199:end); solDT=DT(end-199:end);
         soltime=time((krestart-1)*nnt+1:krestart*nnt); solDT=DT((krestart-1)*nnt+1:krestart*nnt);
         save stabtrNS_end.mat U Udot solDT soltime DT time initdata dateflag t tol viscosity
         fname=sprintf('%s%i','MovieData',krestart);   %  filename for movie data
         save(fname,'U','Udot','soltime');             %  save movie data
         uxmaxs=max(max(abs(U(1:nv,:))));
         uymaxs=max(max(abs(U(nv+1:2*nv,:))));
         uxmins=min(min(U(1:nv,:)));
         uymins=min(min(U(nv+1:2*nv,:)));
      % update extreme value datafile
            load('VelocityScaling.mat');                        %  load previous maximum
            uxmax=max(uxmax,uxmaxs); uymax=max(uymax,uymaxs);   %  new maximum
            uxmin=min(uxmin,uxmins); uymin=min(uymin,uymins);   %  new minimum
            save('VelocityScaling.mat','uxmax','uymax','uxmin','uymin');
         time=longtime; DT=longDT;
         fprintf(' --- CheckPoint')
         if(tout>=1)
            fprintf('\n\n   step  timestep   time     meanKE   vorticity'); 
         end
         k=0;      % reset storage counter  	   
      end 
      n = n+1; k=k+1; 
   else
%%%%%% rejected step     
      nrej = nrej + 1; 
      if(tout==1)
         fprintf(' oops .. step was rejected'); 
      end  
   end
  
%%%%%% compute the next timestep
   dt = dt*(tol/d)^(1/3);
end        %---------------- end of timestep loop

%
%%% finishing touches
fprintf('\nfinished in %4i steps!\n',n-1)
DT = [DT(1:n-1) dt]; time = time(1:n-1); %%%%%%
%% solution is returned for time steps after the final check-point
if(nrej>0) 
   disp([int2str(nrej),' rejections in total: tol = ',num2str(tol)]);
end
% resave the final solution
finalt=t;
initdata=struct('u',u,'ub',ub,'udot',udot,'udotb',udotb,'dt0',dt,...
                'dt',dt*(tol/d)^(1/3),'time',t,'n',n-1,'restart',restart');
fprintf('final time solution data is saved in stabtrNS_end.mat\n');
save stabtrNS_end.mat DT time initdata finalt tol 
if(k>1)
   U = U(:,1:k-1); Udot = Udot(:,1:k-1); 
   uxmaxs=max(max(abs(U(1:nv,:))));
   uymaxs=max(max(abs(U(nv+1:2*nv,:))));
   uxmins=min(min(U(1:nv,:)));
   uymins=min(min(U(nv+1:2*nv,:)));
   krestart=krestart+1;
   [s1,s2]=size(U);
%  soltime=time(end-s2+1:end);
   soltime=time((krestart-1)*nnt+1:end); solDT=DT((krestart-1)*nnt+1:end);
   fname=sprintf('%s%i','MovieData',krestart);    %  filename for movie data
   save(fname,'U','Udot','soltime');              %  save movie data
   % update extreme value datafile
        load('VelocityScaling.mat');
        uxmax=max(uxmax,uxmaxs); uymax=max(uymax,uymaxs);
        uxmin=min(uxmin,uxmins); uymin=min(uymin,uymins);
        save('VelocityScaling.mat','uxmax','uymax','uxmin','uymin');
end
save('MovieSegments','krestart');          %  save number of movie segments
return


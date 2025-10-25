function [DT,U,Udot,time] = stabtrx(A,M,f,uzero,dtzero,tfinal,tol,nstar,info)
%STABTRX  vanilla stabilised TR | extended
%   [DT,U,Udot,time] = stabtrx(A,M,f,uzero,dtzero,tfinal,tol,10,0);
%   input
%          A, M      specified ODE system: M udot + A u = f without BCs
%          uzero     initial condition
%          dtzero    initial timestep
%          tfinal    final time 
%          tol       local accuracy error tolerance   
%          nstar     averaging frequency   
%          info      output switch
%   output
%          DT        timestep history
%          U         solution history
%          Udot      solution time derivative history
%          time      discrete time evolution
%   X-IFISS function: DJS; 15 March 2018
% Copyright (c) 2018 D.J. Silvester, Huining Yang

%--------- setup
nmax=2e5;  % maximum number of timesteps
ub=uzero; nd=length(ub);
T=tfinal; dt0=dtzero;
fprintf('StabTR iteration in %g dimensions ...',nd)
% preallocate array dimensions
DT = zeros(1,50); U = zeros(nd,50); time = zeros(1,50); Udot = zeros(nd,50);

%--------- initialization
ff=f; f=0; %no forcing at t=0
udotb = M\(-A*ub+f);                                % du/dt(0)
ke = sqrt((ub'*(M*ub))); acc = sqrt((udotb'*(M*udotb)));
itntable(nd,0,0,0,ub,udotb,ke,acc)

%--------- first time step
v = (M+ 0.5*dt0*A)\(f + M*udotb - A*ub);               % first TR step
u = ub + 0.5*dt0 *v;
udot = 2*(u-ub)/dt0 - udotb;	                       % du/dt(dt0)
udd = (udot-udotb)/dt0;		                           % second derivative
dt = dt0; t=dt;
r = norm(M*udot+A*u-f);
ke = sqrt((u'*(M*u))); acc = sqrt((udot'*(M*udot)));
itntable(nd,dt0,r,0,u,udot,ke,acc)

n=2;                                                   % set time step index
DT(1:2) = [dt0, dt];
U(:,1) = ub; U(:,2) = u;
Udot(:,1) = udotb; Udot(:,2) = udot;
time(1:2) = [0,dt0];

flag = 0; nrej=0; avflag = 0;  nav=nstar; tstar=tfinal;

%--------- loop until time limit is reached
while t <= T  & flag==0
  if t+dt>T, dt = T-t; flag = 1; end                   % fix final time step
     if n==nmax, flag=1;
     fprintf('\nToo slow -- step limit reached!'), end
  f = (pi*pi*sin(t)+cos(t))*ff;                        % sinusoidal forcing in time
  v = (M+0.5*dt*A)\(f + M*udot - A*u);                 % general TR step
  w = udot + 0.5*dt*udd;                               % AB2 step
  udiff =  0.5*v-w;
  upTR  = u + 0.5*dt*v;
  d = (dt*dt/(3*(dt+dt0)))*sqrt( udiff'*M*udiff);

% timestep control
  if d < 2*((1/0.7)^3)*tol      % time step accepted
    if (t>tstar & avflag==0) | ~rem(n,nav)    %%% smooth by averaging
       dt0 = 0.5*(dt+dt0);
       ub  = 0.5*(u+ub);
       udotb = 0.5*(udot+udotb);
       u = 0.5*(u + upTR);
       udot = 0.5*v;
       t = t + 0.5*dt;
       avflag=1;
       if nav == 1e4, nav = n; end
       if info==1, disp(['Averaging: n = ' int2str(n) ', t = ',num2str(t)]), end
       else                                   %%% take regular step
       dt0 = dt; t = t+dt0;
       ub = u; u = upTR;
       udotb = udot; udot = v - udot;
    end
  udd = (udot-udotb)/dt0;
  r = norm(M*udot+ A*u-f);                      % sanity check
  n=n+1;
    if n > length(time)		                    % allocate more memory
       DT   = [DT   zeros(1,100)];
       U    = [U zeros(nd,100)];   Udot = [Udot zeros(nd,100)];
       time = [time zeros(1,100)];
       end    
  DT(n) = dt; U(:,n) = u;  Udot(:,n) = udot; time(n) = t;
  else   % rejected step
  nrej = nrej + 1; 
  if info==1, disp(['oops .. step ', int2str(n),' rejected']), end
  end
  
% compute the next timestep
dt = dt*(tol/d)^(1/3);
ke = sqrt((u'*(M*u))); acc = sqrt((udot'*(M*udot)));
itntable(nd,t,r,d,u,udot,ke,acc)
end
%--------- end of timestep loop


%--------- finishing touches
fprintf('finished in %3i steps!\n',n)
DT = DT(1:n); U = U(:,1:n); Udot = Udot(:,1:n); time = time(1:n);
if nrej>0, disp([int2str(nrej),' rejections in total: tol = ',num2str(tol)]), end
return
   

function itntable(nd,t,r,d,u,udot,ke,acc)
if nd==1,  % one-dimensional problem
   if t==0;
   fprintf('\n t          u           ke        d\n'),end
   fprintf(' %6.2e  %9.4e  %9.4e  %9.4e\n',t,u(1),sqrt(u'*u),d)
end
if nd>1,  % nd-dimensional problem
if t==0; fprintf('\n     t        ke         acc          d          res  \n'), end
fprintf(' %7.3e  %9.4e  %9.4e   %9.4e  %9.4e\n',t,ke,acc,d,r)
end
return

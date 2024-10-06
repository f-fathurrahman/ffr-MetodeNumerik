% PURPOSE: Compute and plot the time evolution 
%  of a Korteweg-De-Vries (KdeV) soliton, to see how
%  in nonlinear problems, solutions at the
%  end may have higher frequency components
%  than at the beginning
% USAGE: It is possible to reduce the spatial
%  resolution ds, then npoints should be 
%  increased to retain the overall shape:
%     ds=ds/2 ->  npoints=npoints*2
%  IF ds is reduced, dt has to be decreased
%  at least by the same amount, or stronger
%  (vonNeumann stability):
%     ds=ds/2 -> dt=dt/2 
%     ds=ds/3 -> dt=dt/3
%     ds=ds/4 -> dt=dt/5 etc.
% ALGORITHM: Finite difference approximation
%  for the KdeV-equation:
%
%  1. Approximate
%     du            du        d^3 u
%    ----- + eps u ---- + mu -------  = 0
%     dt            dx        dx^3
%  approximate the first time derivative df/dt 
%  with centered differences in second order as
%     du(i,j)     u(i,j+1)-u(i,j-1) 
%    --------- = -------------------  
%       dt                dt
%  ("three-point formula", as the points j,j-1,j+1 
%  are involved) likewise the spatial derivative
%     du(i,j)    u(i+1,j)-u(i-1,j) 
%    --------- = -------------------  
%        dx              dx
%
%  2. Approximate the third derivative with the 
%   "five point formula" (as the points 
%   j,j-1,j+1,j-2,j+2 are involved)
%     d^3u(i,j)   u(i+2,j)-2u(i+1,j)+2u(i-1,j)-u(i-1,j) 
%    --------- = ------------------------------------------  
%        dx^3                    2 dx^3
%    
%  3. Rewrite the terms so that an explicit integrator
%    u(i,j+1)=u(i,j) + ..... results
%
%
% CAVEAT: 
%  - Small "wriggles" develop near the boundaries and 
%   travel "inside", because the discretization at the 
%   boundary is not so accurate/stable
%  - Because the problem is non-linear, the initial 
%   condition cannot be modified arbitrarily without 
%   destroying the stability of the whole solution
% REVISION HISTORY:
%  The computation of the third derivatives  on 
%  the left and right boundary has been modified
%  compared to the book's listing to reduce the noise;
%  The same formula as "inside" the system is used,
%  while the value to "the left" and to "the right" of
%  the boundary is extrapolated.
% TODO: Instead of using extrapolated values in the
%  third order finite difference "on the inside" for the
%  boundary, it would be more satisfying to use higher-order
%  assymmetric formulae which do not need extrapolated
%  values
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis


clear all,  format compact
disp('Solver the Korteweg-DeVries Equation')
ntime=9000 % number of timesteps dt
npoints=131 % number of intervals of length ds
dt=0.025 % size of the timestep 
mu=0.1  % Prefactor for the term with the third derivative
eps=0.2 % Prefactor for the term with the nonlinearity
ds=0.4  % Grid-size
u(:,1)=0.5*(1-tanh(.2*ds*([1:npoints]-1)-5)); % Initial state

if max(u)>1.1
  warning('!!!!!!!!!!!!!!!!!')
  disp('Your initial values are larger than in the original')  
  disp('version. Dont forget to rescale the graphics axis, ')  
  disp('else part of the solution cannot be seen!!')  
end    
if (ds<0.02)
  warning('!!!!!!!!!!!!!!!!!')
  disp('When you reduce the spatial resolution ds to,')
  disp('ds/n, you should reduce dt to dt/(n+....)')
end    
    
% Left Boundary for the next timesteps:
u(1,2)=u(1,1);  
u(1,3)=u(1,1);
% Right Boudnary for the next timesteps:
u(end,2)=u(end,1);
u(end,3)=u(end,1);


fac=mu*dt/(ds^3.0);
time=dt;
for i=2:npoints-1 % First timestep@
  a1=eps*dt*(u(i+1,1)+u(i,1)+u(i-1,1))/(ds*6);
  if ((i>2)&(i<=npoints-2)) % Inside the system
    a2=u(i+2,1)+2*u(i-1,1)-2*u(i+1,1)-u(i-2,1);
  elseif ((i==2)|(i==npoints-1)) % On the left boundary
    a2=u(i-1,2)-u(i+1,1);
  end
  a3=u(i+1)-u(i-1,1);
  u(i,2)=u(i,1)-a1*a3-fac*a2/3;
end

for j=1:ntime % time evolution
  for i=2:npoints-1
    a1=eps*dt*(u(i+1,2)+u(i,2)+u(i-1,2))/(ds*3);
    if ((i>2)&(i<=npoints-2)) % Inside the system
      a2=u(i+2,2)+2*u(i-1,2)-2*u(i+1,2)-u(i-2,2);
    elseif (i==2) % On the left boundary
% Formulae in the book's listing            
%      a2=u(i-1,2)-u(i+1,2) is replaced with a more 
% stable version where the value left from the boundary
% is extrapolated
      uim22=u(i-1,2)+ds*(u(i-1,2)-u(i,2));
% and used in the finite differenc approximation for
% "inside the system" 
      a2=u(i+2,2)+2*u(i-1,2)-2*u(i+1,2)-uim22;
    elseif (i==npoints-1) % On the right boundary
% Formulae in the book's listing    
%      a2=u(i-1,2)-u(i+1,2) is replaced with more 
% stable version where the value left from the boundary
% is extrapolated
      uim2p=u(i+1,2)-ds*(u(i+1,2)-u(i,2));
% and used in the finite differenc approximation for
% "inside the system"
      a2=u(i+1,2)+2*u(i-1,2)-2*u(i+1,2)-u(i-2,2);
    end
    a3=u(i+1,2)-u(i-1,2);
    u(i,3)=u(i,1)-a1*a3-2*fac*a2/3;  
  end
  u(:,1)=u(:,2); % work only with the previous
  u(:,2)=u(:,3); % and the current time-step

  if (20*round(j/20)==j)
     clf 
    plot(ds*[1:npoints],u(:,1),'k-')
    xlabel('x')
    ylabel('Amplitude')
    axis([1 ds*npoints -.5 2])
    drawnow
  end
end

return

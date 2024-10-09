%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2023                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Initial condition for the vortex dipole       %
% computes the velocity field for a single vortex %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,v]=NSE_F_init_vortex(Lx,Ly,x,y,xv,yv,lv,uin,vin,pm)

global nxm nym dx dy
global im ip jp jm ic jc

% xv, yv  : coordinates of the vortex center
% uin,vin : initial translation velocity

            % radius of the vortex
%
            % intensity $\psi_0$
psi0=0.1;
            % stream-function and velocity field
for jy=1:nym
for jx=1:nxm
  uloc=(x(jx)-xv)^2+(y(jy)-yv)^2;
  uloc=pm*psi0*exp(-uloc/lv^2);
  u(jx,jy)=uin-2.d0*(y(jy)-yv)*uloc/lv^2;
  v(jx,jy)=vin+2.d0*(x(jx)-xv)*uloc/lv^2;
end
end

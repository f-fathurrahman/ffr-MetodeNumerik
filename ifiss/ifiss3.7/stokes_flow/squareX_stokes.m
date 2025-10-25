%SQUAREX_STOKES set up two-field approximation in unit square domain
%   IFISS scriptfile: DJS; 4 August 2023
% Copyright (c) 2023 D.J. Silvester, J. Pestana
clear variables
%% define geometry
pde=3; domain=1; enclosed=1;
cavity_domain
load cavity_grid.mat
%
%% set up matrices
qmethod=5;
[x,y,xy,xyp,mp,map] = q2q1gridx(x,y,xy,mv,bound);
[A,B1,B0,Q1,Q0,G,Bx,By,f,g1,g0] = stokes_q2q1q0(xy,xyp,mv,mp);
%
% construction of pressure mixing matrix R
rlocal=diag(Q0)*[1,1,1,1]; np0=length(g0); np1=length(g1);
R=sparse(np0,np1);
% sparse assembly ...
   for pp=1:np0
      R(pp,mp(pp,:))=Q0(pp,pp)*[1,1,1,1]/4;
   end
Q=[Q1,R';R,Q0];
fprintf('Augmentation of Taylor-Hood approximation completed\n')

gohome
cd datafiles
save squareX_stokes_nobc.mat ...
pde domain enclosed qmethod grid_type A B1 B0 G Q1 Q0 Q Bx By f g1 g0
save squareX_stokes_nobc.mat mv mp map xy xyp mbound bound x y  -append

fprintf('system matrices saved in squareX_stokes_nobc.mat ...\n')

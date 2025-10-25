%NEWCHANNEL_DOMAIN rectangular shaped domain Q2 grid generator
%   newchannel_domain;
% 
% Prior to calling the script file the variable outbnd may
% be set to the outflow length.
% natural ouflow boundary condition  is applied on right-hand edge
% grid defining data is saved to the file: channel_grid.mat
%   IFISS function: HCE; 24 December 2014, DJS; 9 January 2019
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified (c) 2023 M.D. Mihajlovic

%clear variables
fprintf('\n\nGrid generation for extended channel domain.\n')
if exist('outbnd','var')==0,
fprintf('setting channel length to default\n'), outbnd=1; end
nc=default('grid parameter? (default is 4)',4); nc=nc-1; %%% compatibility
if nc<1, error('illegal parameter choice, try again.'), end
grid_type=1;
ny=2^nc; npy=ny/2; nqy=ny/4;                   %%% modified from channel_domain
nx=ny*(outbnd/2);  npx=nx/2;  nqx=nx/4;        %%% to handle rectangular domain
nx = min(ny*(outbnd/2),ny*100/2);   npx = nx/2;   nqx = nx/4;
%
%
%% compute (x,y) coordinates of vertices
yy=(1/npy:1/npy:1);
ypos=[0,yy]; yneg=-yy(length(yy):-1:1);
yyy=[yneg,ypos];
xxx=(-1:(outbnd/2)/npx:outbnd);
%
[xy,mv,left,right,bottom,top,mbound,x,y] = grid_xblock(xxx,yyy);
outbc=0;
% exclude outflow boundary data
bound=unique([left;top;bottom]);
kk=find(mbound(:,2)~=2); mbound=mbound(kk,:);
%
%% specify boundary information for graphics
% bndxy: (x,y)-coordinates of vertices that define the domain and
%         obstacle(s) 
% bnde: boundary edges (node1 node2 1(for dirichlet)/0(for neumann))
% obs: obstacles (node1 node2 node3 node4)
% sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
% 'obs' and/or 'sbnde' can be absent if there is no obstacle in the problem
% and/or only uniform grid is needed
bndxy = [-1,-1; outbnd,-1; outbnd,1; -1,1];
bnde = [1,2,1; 2,3,0; 3,4,1; 4,1,1];
obs = [];
sbnde = [];

%%
gohome
cd datafiles
save channel_grid.mat mv xy bound mbound outbnd grid_type outbc x y bndxy bnde obs 
fprintf('generated datatfile: channel_grid.mat \n')
%
%% Plot the domain geometry for movies
%
figure(2)
xp=[bndxy(1,1),bndxy(2,1)]; yp=[bndxy(1,2),bndxy(2,2)];   %  bottom horizontal boundary
plot(xp,yp,'-k','Linewidth',3); hold on
xp=[bndxy(2,1),bndxy(3,1)]; yp=[bndxy(2,2),bndxy(3,2)];   %  right vertical boundary
plot(xp,yp,'-g','Linewidth',3); hold on
xp=[bndxy(3,1),bndxy(4,1)]; yp=[bndxy(3,2),bndxy(4,2)];   %  top horizontal boundary
plot(xp,yp,'-k','Linewidth',3); hold on
xp=[bndxy(4,1),bndxy(1,1)]; yp=[bndxy(4,2),bndxy(1,2)];   %  left vertical boundary
plot(xp,yp,'-r','Linewidth',3); hold on
axis([-1.2 outbnd+0.2 -1.2 1.2]);
axis equal; axis off;
channel_dom=gca;
save channel_plot.mat channel_dom
close(2)
return

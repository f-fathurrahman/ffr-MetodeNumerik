close all;
clear variables;

global dx dy Lx Ly;
global nxm nym;
global ip im jp jm ic jc;

Lx = 2.0;
Ly = 1.0;
nx = 5;
ny = 5;

nxm = nx - 1;
nym = ny - 1;  %number of cells
dx = Lx/nxm;
dy = Ly/nym;

ic = 1:nxm;
jc = 1:nym;  % indices of cells

xc = (ic-1)*dx;
yc = (jc-1)*dy;% primary grid

xm = (ic-0.5)*dx;
ym = (jc-0.5)*dy;  % mid-cell grid

ip = ic + 1;
ip(nxm) = 1;

jp = jc + 1;
jp(nym) = 1;  % indices for periodicity

im = ic - 1;
im(1) = nxm;

jm = jc - 1;
jm(1) = nym;  % cell 1 = cell nxm+1
                             
[xx, yy] = meshgrid(xm, ym);
% why transpose? to make the shape to be (nxm,nym)
xx = xx';
yy = yy'; % centers of the cells for visualization

% initial velocity field (icas=1)
u  = NSE_F_init_KH(Lx, Ly, xc, ym, 1, Ly/4, 20, 0.25, 0.5*Lx);
sca = NSE_F_init_KH(Lx, Ly, xm, ym, 1, Ly/4, 20, 0.00, 0.5*Lx);

% Purpose: Driver for DEMround2Dnorot for the 
%  two-dimensional simulation of circular particles
% CAVEAT: 
%  - There are no rotational degrees
%    of freedom for the particles
%  - The visualization may not be smooth
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

clear all
format compact
clf

n_part=18
% initialize radius and mass
global rad, rad(1:n_part)=0.5;
global m, m(1:n_part)=1;
global E, E=10000; % Young modulus
global lmaxx, lmaxx=n_part/2+1;
global lminx, lminx=0;
global lmaxy, lmaxy=n_part/2+1;
global lminy, lminy=0;
global g, g=-9.81;

% initialize positions and velocities=0
rand('seed',5)
r0_x=2*mod([1:n_part],5)+1;
v0_x=rand(size(r0_x))-0.5;
r0_y=sort(r0_x);
v0_y=rand(size(r0_y))-0.5;
y0(1:4:4*n_part-3)=r0_x;
y0(2:4:4*n_part-2)=v0_x;;
y0(3:4:4*n_part-1)=r0_y;
y0(4:4:4*n_part)=v0_y;;

t_end=5;
[t,y]=ode113('DEMround2Dnorot',[0:0.005:t_end],y0);

figure(1)
[X,Y,Z]=cylinder(rad,55); % outline for plotting

hold on
for i=1:length(t)
  clf
  hold on
  for i_part=1:n_part
    c1=X(i_part,:)+y(i,4*i_part-3);
    c2=Y(i_part,:)+y(i,4*i_part-1);
    fill(c1,c2,[.3 .3 .3]) % plot particles
  end  
  a=quiver(y(i,1:4:4*n_part),y(i,3:4:4*n_part),...
         y(i,2:4:4*n_part),y(i,4:4:4*n_part));
  set(a,'Color',[0 0 0]) % plot velocity with black arrows
  axis equal
  axis([lminx-.5 lmaxx+.5 lminy-.5 lmaxy+.5])
  drawnow  
end    

return
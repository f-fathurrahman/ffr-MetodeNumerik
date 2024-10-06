function [polygonX,polygonY,initialX,initialY,initialF,...
  polygonSide,polygonisfree,m,mi]=demInitial
% PURPOSE: 
% 1. Read in the geometry of the walls
% 2. Initialize the particles
% 3. Compute the masses, moments of inertia 
%    of the particles and wall
% 4. Concatenate the data for the respective particles 
%    and walls so that both can be treated together 
% OUTPUT: 
%  polygonX,polygonY contains the corner positions,
%  initialX,initialY,initialF contains the positions
%  and their higher time derivatives 
% CAVEAT: 
%  - If two particles are initialized with overlaping
%   bounding boxes, the neighborhood algorithm may fail 
%   to detect this pair of interacting particles
%  - Wall particles are considered to be "infinitely
%   heavy and are rescaled to 1000 times the 
%   largest (physical) mass of any particle or wall
% PERFORMANCE CONSIDERATION: The MATLAB-editor gives
%  warnings for several variables that the size
%  changes in every loop; because this happens only
%  for the initialization, no performance optimization
%  was attempted
% TODO: Write an additional routine which makes sure
%  that when a new particle is initialized, it has no 
%  overlap with the previous ones; This can be done by
%  the call of the overlap-routine for the newly initialized
%  particle and a loop over all previously initizalied 
%  particles 
% LITERATURE: sec. 7.4
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014
% 

global rhoPoly
global rhoWall
global gravity

% Read in walls first: From the position of the
% walls, that domain can be computed in which the 
% particles should be initialized 
walldatain=fopen('walldata.dat')
dummy=textscan(walldatain,'commentStyle','%');
dummy=textscan(walldatain,'commentStyle','%');
wallNum=textscan(walldatain,'%d1','commentStyle','%');
wallNum=cell2mat(wallNum);
for i=1:wallNum
  numcorn=textscan(walldatain,'%d1','commentStyle','%');
  numcorn=cell2mat(numcorn)
  wallSide(i)=numcorn;
  isfree=textscan(walldatain,'%d1','commentStyle','%');
  isfree=double(cell2mat(isfree))
  wallisfree(i,1)=isfree;
  for j=1:wallSide(i)
    wallpos=fscanf(walldatain,'%g %g\n',[1 2])
    wallX(j,i)=wallpos(1)
    wallY(j,i)=wallpos(2)  
  end
end  


% Initialize the maximal particle radius with 
% the width of half the left wall
particle_radius=.5*max(wallX(:,2))-min(wallX(:,2));

% for walldata.dat, 
% Wall 1 is the floor,
% Wall 2 is the left, Wall 3 the right wall,
% Wall 4 is the ceiling,
% Wall 5 and 6 are the left and right slopes 
% of the hopper.
% Initialize the centers of mass of the
% particles so that they are in the
% upper part of the hopper:
xmin=max(wallX(:,2))+1.7*particle_radius;
xmax=min(wallX(:,3))-1.7*particle_radius;
ymin=max(wallY(:,5))+1.7*particle_radius;
ymax=min(wallY(:,4))-1.7*particle_radius;
nx=floor((xmax-xmin)/ particle_radius/2); 
ny=floor((ymax-ymin)/ particle_radius/2); 

xpos=linspace(xmin,xmax,nx);
ypos=linspace(ymin,ymax,ny);

particleNum=0;
for ix=1:nx
  for iy=1:ny
    particleNum=particleNum+1;
% number of corners between 4 and 8    
    ncorner=4+round(rand*3+.5);
% angles for the corners
    angles=2*pi/ncorner*[0:ncorner-1];
% Elongated particles, ellipse radius:
    ar=(.8+.2*rand)*particle_radius;
    br=(.8+.2*rand)*particle_radius;
% Center of mass 
    for k=1:ncorner
      particleX(k,particleNum)=xpos(ix)+cos(angles(k))*ar;
      particleY(k,particleNum)=ypos(iy)+sin(angles(k))*br; 
    end
    particleSide(particleNum)=ncorner;
    particleisfree(particleNum,1)=1;
    % Initialize velocities as 0
    VX(particleNum,1)=0;
    VY(particleNum,1)=0;
    VF(particleNum,1)=0;
  end
end

% Initialize masses, mom. of inertia, centers
% For Walls
[mwall,miwall,wallCOMF,wallCOMX,wallCOMY]=...
 mass_momentinertia(wallSide,wallX,wallY,rhoWall);
% For Particles
[m,mi,particleCOMF,particleCOMX,particleCOMY]=...
 mass_momentinertia(particleSide,...
 particleX,particleY,rhoPoly);
% Because the walls should be "infinitely heavy",
% their masses are set to "large" values (1000 times 
% larger than the heaviest particle/wall)
maxm=max([m mwall]);
maxmi=max([mi miwall]);
for i=1:wallNum
  m(i) =maxm;
  mi(i)=maxmi;
end

% Combine data structures for particles and walls
maxSide = max([particleSide wallSide]); 
polygonNum = particleNum+wallNum;
m=[m mwall];
mi=[mi miwall];
polygonX=zeros(maxSide,polygonNum);
polygonY=zeros(maxSide,polygonNum);
polygonSide=[wallSide particleSide];
polygonisfree=[wallisfree; particleisfree];

% Copy walls into the particle structure
for i=1:wallNum
  tempSide = wallSide(i);
  polygonX(1:tempSide,i) = wallX(1:tempSide,i);
  polygonY(1:tempSide,i) = wallY(1:tempSide,i);
end
iparticle=0;
% Copy particles into the particle structure
for i=wallNum+1:wallNum+particleNum
  iparticle=iparticle+1;
  tempSide=particleSide(iparticle);
  polygonX(1:tempSide,i)=particleX(1:tempSide,iparticle);
  polygonY(1:tempSide,i)=particleY(1:tempSide,iparticle);
end


% Variables for predictor the predictor corrector:
% up to six values, the positions and its 5 time
% derivatives:
% 1. Center of mass
% 2. Velocity
% 3. Acceleration
% 4. Jerk
% 5. First derivitive of jerk.
% 6. Second derivitive of jerk.
initialX = zeros(particleNum,6);
initialY = zeros(particleNum,6);
initialF = zeros(particleNum,6);
for i=1:wallNum
% wall velocities and their time derivatives are 0,
% if walls moving at constant velocity
% are desired, changes must also be made
% in the predictor and the corrector, where 
% in the current version the wall velocities are 
% set to zero.
  initialX(i,:)=[wallCOMX(i) 0 0 0 0 0]; % x-position
  initialY(i,:)=[wallCOMY(i) 0 0 0 0 0]; % y-position
  initialF(i,:)=[wallCOMF(i) 0 0 0 0 0]; % orientation
end
for i=1:particleNum
  initialX(wallNum+i,:)=[particleCOMX(i) VX(i) 0 0 0 0];        
  initialY(wallNum+i,:)=[particleCOMY(i) VY(i) -gravity 0 0 0]; 
  initialF(wallNum+i,:)=[particleCOMF(i) VF(i) 0 0 0 0];        
end

return
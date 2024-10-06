% PURPOSE: Companion program for 
%  Understanding the Discrete Element Method: 
%  Simulation of Non-Spherical Particles for 
%  Granular and Multi-body Systems
%  Hans-Georg Matuttis, Jian Chen
%  Wiley 2014
%  to which the page- and equation numbers refer
% THIS PROGRAM AND ITS FUNCTIONS COME WITHOUT WARRANTY!
% USAGE: Driver for the discrete element method
%  HGM_DEM_2D with polygonal particles
%  Simulation of the outflow through a hopper 
% CAVEATS: The purpose of this program is not to
%  supply a perfectly running code (though it should
%  be sufficiently bug-free to use it for simulations),
%  but to illustrate the concepts from the book, including
%  various programming approaches, documentation examples.
% AUTHORS: Shi Han Ng, Hans-Georg Matuttis
% TODO:
%  The items in the TODO-list are proposals how the respective
%  functions could be improved 
%   
%  Eliminate particle number in the argument list
%  Input page- and formula-numbers
%  Input the final date of the revision history
%  when all functions are satisfyingly documented

clear all
clear neighborhood_algorithm % Clear persistent variables in neighborhood_algorithm.
format compact

startingTime = tic; % cputimestart = cputime;

global gravity rhoPoly rhoWall emodul mu gamma  

% set seed to obtain the same configuration 
% in every program run: only for debugging
rand('seed',3) 

% Time-variables
t     = 0.0            % Initial simulation time
tmax  = .20;           % Final simulation time 
dt    = 1e-5;          % Size of the timestep 
nt=round((tmax-t)/dt); % Number of steps 
everyn= 100;           % Number of steps between graphical output

% DEM-parameters
gravity   = 9.81;    % [m/s^2] Standard acceleration due to free fall.
rhoPoly   = 5000.0;  % [kg/m^3] Density of particle.
rhoWall   = rhoPoly; % [kg/m^3] Density of wall.
emodul    = 1.0e+6;  % Young modulus. 
mu        = 0.3;     % Coefficient of friction.
gamma     = 1.5;     % Damping constant.

% Initialization of the particles; initial coordinates named 
% "corrected" to avoid recopying with other variable names
% for the main loop
[particleX,particleY,correctedX,correctedY,correctedF,...
 particleSide,freeParticle,m,mi]=demInitial;

% Maximum stepsize according to the natural frequency of 
% the lightest particle, p.240.
% Analogous to the discussion on 
% p. 291 for the three-dimensional case
MAXDT = 0.1*sqrt( min(m)/emodul )*pi;
if ( dt > MAXDT )
  error('The stepsize is too large.'); 
elseif(dt<5*MAXDT)
  disp(['Warning! The timestep of' num2str(dt) ])  
  disp('seems to be considerably smaller than the')
  disp(['maximally allowed timestep of' num2str(MAXDT) ])  
  dummy=input('contiue? Press any key') ;
end

% Graphical output of the initial configuration,
% to make sure that particles are where they should be
figure(1), clf, hold on
draw_particle(particleSide,particleX,particleY);
axis image

% Main loop
disp('Starting time integration')
for it=1:nt
% Increment the simulation time.
  t=t+dt;
% Graphics output every everyn steps
  if (everyn*round(it/everyn)==it)    
    it
    figure(1), clf, hold on
    draw_particle(particleSide,particleX,particleY)
    axis image
  end

% Predictor step 
  [PX,VX,PY,VY,PF,VF] = gearPredictor(dt,freeParticle,...
  correctedX,correctedY,correctedF);

% Update the vertices according to the predictor coordinates
  [particleX,particleY]=updatevertices(particleSide,...
  correctedX(:,1),correctedY(:,1),correctedF(:,1),...
  PX,PY,PF,particleX,particleY);

% Neighborhood algorithm
% (initialization inside the function)
  [col_list,col_list_len]=neighborhood_algorithm(...
  particleSide,particleX,particleY);

% Interaction computation 
  [acx,acy,acf]=deminteraction(col_list,col_list_len,...
  dt,particleX,particleY,particleSide,...
  PX,PY,VX,VY,VF,freeParticle,m,mi);   

% Corrector-step of the Gear Predictor-Corrector 
  [correctedX,correctedY,correctedF]=...
  gearCorrector(acx,acy,acf,dt,freeParticle);
    
end % end of time evolution loop

disp('Finished without error!')

return
clear all; clc; close all; clear workspace
input = init();
global nDIM;

if nDIM == 1
   % Compute scattered wave field at a receiver above and below the slab
     disp('Running WavefieldSctSlab');    WavefieldSctSlab;
     
elseif nDIM == 2
   % Compute scattered wave at receivers around the circular cylinder
     disp('Running WavefieldSctCircle');  WavefieldSctCircle;  
    
elseif nDIM == 3
   % Compute scattered wave at a number of receivers around the sphere
     disp('Running WavefieldSctSphere');  WavefieldSctSphere;  
     
end % if

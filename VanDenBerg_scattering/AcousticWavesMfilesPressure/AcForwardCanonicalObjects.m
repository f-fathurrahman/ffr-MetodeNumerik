clear all; clc; close all; clear workspace
input = initAC();
global nDIM;

if nDIM == 1; 
  % Compute scattered acoustic field at receivers above/below the slab
  disp('Running AcousticSctWavefieldSlab');    AcousticSctWavefieldSlab;
     
elseif nDIM == 2; 
   % Compute scattered acoustic field at receivers around circular cylinder
   disp('Running AcousticSctWavefieldCircle'); AcousticSctWavefieldCircle;
    
elseif nDIM == 3; 
   % Compute scattered acoustic field at receivers around sphere
   disp('Running AcousticSctWavefieldSphere'); AcousticSctWavefieldSphere;

end % if

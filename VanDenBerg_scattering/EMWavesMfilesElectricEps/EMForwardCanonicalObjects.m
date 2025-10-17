clear all; clc; close all; clear workspace
input = initEM();
global nDIM;

if nDIM == 2 
   % Compute scattered acoustic field at receivers around circular cylinder
   disp('Running EMsctCircle'); EMsctCircle;
    
elseif nDIM == 3 
   % Compute scattered acoustic field at receivers around sphere
   disp('Running EMsctSphere'); EMsctSphere;

end % if
 
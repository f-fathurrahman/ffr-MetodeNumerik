clear all; close all; 
%   SYNTAX 
%   example_decimation_reducepatch
%   DESCRIPTION 
%   This script implements edge decimation using built-in MATLAB function
%   REDUCEPATCH. Non-manifold meshes are frequently produced
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open'); load(FileName, '-mat');
[t P] = reducepatch(t, P, 700); % reduce to 500 or so triangles
[P t] = fixmesh(P, t);

X = reshape(P(t', 1),[3, size(t, 1)]);
Y = reshape(P(t', 2),[3, size(t, 1)]);
Z = reshape(P(t', 3),[3, size(t, 1)]);
patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);
light('Position', [0 1 0], 'Style', 'local'); axis 'equal';  axis 'tight'; view(0, 0); grid on;
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');
    
NumberOfTirangles = size(t, 1)
Q = min(simpqual(P, t))
edges = meshconnee(t);
temp = P(edges(:, 1), :) - P(edges(:, 2), :);
minedgelength = min(sqrt(dot(temp, temp, 2)));

NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
save(NewName, 'P', 't');

[t, flag1]  = checkmanifold(t);
flag2       = checkintersection(P, t);
if flag1 error('The mesh is not manifold'); end
if flag2 error('The mesh has self-intersections'); end

    






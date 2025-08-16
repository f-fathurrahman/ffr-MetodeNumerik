clear all; close all; 
%   SYNTAX 
%   example_decimation_reducepatch_soft
%   DESCRIPTION 
%   This script implements edge decimation using built-in MATLAB function
%   REDUCEPATCH. An attempt is made to vary the number of triangles in order
%   to create a manifold mesh. Usually does not work
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open'); load(FileName, '-mat');
P0 = P; t0 = t;
for m = 500:100:2000
    [t P] = reducepatch(t0, P0, m); 
    [P t] = fixmesh(P, t);
    
    NumberOfTirangles = size(t, 1)
    Q = min(simpqual(P, t))
    edges = meshconnee(t);
    temp = P(edges(:, 1), :) - P(edges(:, 2), :);
    minedgelength = min(sqrt(dot(temp, temp, 2)));
 
    [t, flag1]  = checkmanifold(t);
    flag2       = checkintersection(P, t);
    if flag1|flag2 
        disp('Mesh check failed') 
    else
        X = reshape(P(t', 1),[3, size(t, 1)]);
        Y = reshape(P(t', 2),[3, size(t, 1)]);
        Z = reshape(P(t', 3),[3, size(t, 1)]);
        patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);
        light('Position', [0 1 0], 'Style', 'local'); axis 'equal';  axis 'tight'; view(0, 0); grid on;
        xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');
        NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
        save(NewName, 'P', 't');
        break;
    end
end

    






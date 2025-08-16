%   viewer_combo
%   DESCRIPTION
%   This script displays two meshes from two *.mat P-t files and major mesh
%   parameters: the number of triangles, minimum triangle quality, and
%   minimum edge length simultaneously
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open'); load(FileName, '-mat');  
edges  = meshconnee(t); 
temp = P(edges(:, 1), :) - P(edges(:, 2), :);
minedgelength1 = min(sqrt(dot(temp, temp, 2)))
Q1 = min(simpqual(P, t))
P1 = P; t1 = t;
    
FileName = uigetfile('*.mat','Select the tissue mesh file to open'); load(FileName, '-mat');
edges  = meshconnee(t); 
temp = P(edges(:, 1), :) - P(edges(:, 2), :);
minedgelength2 = min(sqrt(dot(temp, temp, 2)))
Q2 = min(simpqual(P, t))
P2 = P; t2 = t;

patch('vertices', P1, 'faces', t1, 'EdgeColor', 'k', 'FaceAlpha', 1.0,'FaceColor', 'c');
patch('vertices', P2, 'faces', t2, 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', 'y');
axis 'equal';  axis 'tight'; view(0, 0); grid on;
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');
title('Enclosed nodes: mesh #1 - cyan; mesh #2 - yellow');
    





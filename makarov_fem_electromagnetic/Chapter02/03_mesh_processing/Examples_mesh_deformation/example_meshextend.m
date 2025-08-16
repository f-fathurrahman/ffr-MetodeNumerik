clear all; close all;
%   SYNTAX 
%   example_meshextend
%   DESCRIPTION 
%   This script performs mesh extension/shrinkage in the direction of the
%   normal vector. May be applied to selected nodes
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the mesh file');
load(FileName);

P0 = P; t0 = t;
%  normal direction in mm (positive for extension, negative for shrinkage)
alpha = 0.5; 
%   Find outer normal vectors
normals = meshnormals(P, t);
%   Find triangles attached to every vertex (cell array)
si = meshconnvt(t);
%   Move selected nodes in the direction of the inner normal vector
%   taking into account the local topology (important)
index = 1:length(P);
Ptemp = P(index, :);
for m = 1:length(index)
    n = index(m);
    averagenormal = sum(normals(si{n}, :), 1)/length(si{n});
    norm          = sqrt(dot(averagenormal, averagenormal));
    tangent       = norm;  
    Ptemp(m, :) = Ptemp(m, :) + alpha*(averagenormal/norm)/tangent^1.0;
end 
P(index, :) = Ptemp;    

if alpha>0
    %   View the original mesh
    patch('vertices', P0, 'faces', t0, 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', 'y');
    %   View the deformed mesh
    patch('vertices', P, 'faces', t, 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', 'c');
else
    %   View the original mesh
    patch('vertices', P0, 'faces', t0, 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', 'c');
    %   View the deformed mesh
    patch('vertices', P, 'faces', t, 'EdgeColor', 'k', 'FaceAlpha', 0.5,'FaceColor', 'y');
end

axis 'equal';  axis 'tight'; grid on
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm'); view(90, 0)

normals = meshnormals(P, t);
save(strcat(FileName(1:end-4),'_5mod'), 'P', 't', 'normals');

[t, flag1]  = checkmanifold(t);
flag2       = checkintersection(P, t);
if flag1 error('The mesh is not manifold'); end
if flag2 error('The mesh has self-intersections'); end








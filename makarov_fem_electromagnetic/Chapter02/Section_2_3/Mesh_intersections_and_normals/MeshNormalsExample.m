clear all; close all;
%   SYNTAX
%   MeshNormalsExample
%   DESCRIPTION
%   This script computes and displays the outer normal unit vectors for any
%   2 manifold mesh
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%%  Input structure (*.mat file)
[FileName, PathName] = uigetfile('*.mat','Select the file for normals');
load(FileName);

normals = MeshNormals(P, t);

Center = (P(t(:, 1),:)+P(t(:, 2),:)+P(t(:, 3),:))/3;

%% Plot all patches of the first mesh
Patch = patch('Faces', t, 'Vertices', P, 'EdgeColor', 'r', 'FaceAlpha',1.0,'FaceColor','g');
hold on;
axis 'equal';
xlabel('x, m'); ylabel('y, m'); view(30,60); grid on;
quiver3(Center(:,1), Center(:,2), Center(:,3), normals(:,1), normals(:,2), normals(:,3), 1);

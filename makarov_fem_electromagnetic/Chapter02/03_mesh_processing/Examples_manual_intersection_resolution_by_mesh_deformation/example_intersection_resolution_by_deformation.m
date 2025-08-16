clear all; close all; clc; 
warning off; %  to avoid globals
%   SYNTAX
%   example_intersection_resolution_by_deformation.m
%   DESCRIPTION
%   This scrip deforms meshe(s) by manually moving nodes with given indexes
%   to avoid intersections (nodal check only)

%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

scrsz = get(0,'ScreenSize');

%%  Define box size for a detailed view (important; try to reduce)
delta    = 20;  % in mm
% delta    = 10;  % in mm


%%  Mesh A (*.mat file)
FileName = uigetfile('*.mat','Select the master mesh file to open');
load(FileName, '-mat');
P1 = P; t1 = t;
str.MasterMeshSize = size(t1, 1);
str.MasterMeshQ = min(simpqual(P1, t1));

%% Custom operations with nodes (modify P1) 
% index = [49 70 86 44 64 84 99 47 65 85];
% P1(index, 2) = P1(index, 2) + 10; 
% index = [33 57 83 100 104 23 44 63 80 41];
% P1(index, 2) = P1(index, 2) + 5; 

%% Save the modified file
% NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
% save(NewName, 'P', 't');

%%  Mesh B (*.mat file)
FileName = uigetfile('*.mat',strcat('Select the slave mesh file to open ',FileName));
load(FileName, '-mat');
P2 = P; t2 = t;
str.SlaveMeshSize  = size(t2, 1);
str.SlaveMeshQ = min(simpqual(P2, t2));

%%  Find all nodes of mesh #2 within mesh #1
[in2, trinumber1] = meshinoutstatus(P1, t1, P2);
in2 = find(in2); 

%%  Find all nodes of mesh #1 not within mesh #2
[in1, trinumber2] = meshinoutstatus(P2, t2, P1);
in1 = find(in1); 

if isempty(in1)&isempty(in2)
    a = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
    patch('Faces', t1, 'Vertices', P1, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c');
    patch('Faces', t2, 'Vertices', P2, 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'FaceColor', 'y');
    axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
    xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid off;
    display('Both meshes do not have nodes inside each other (but they may still intersect!)')
    return;
end

%%  Plot original meshes and all enclosed nodes (select the box to zoom in)
center1  = 1/3*(P1(t1(:, 1), :) + P1(t1(:, 2), :) + P1(t1(:, 3), :));
center2  = 1/3*(P2(t2(:, 1), :) + P2(t2(:, 2), :) + P2(t2(:, 3), :));
xbox1     = [min(P1(in1, 1)) max(P1(in1, 1))];
ybox1     = [min(P1(in1, 2)) max(P1(in1, 2))];
zbox1     = [min(P1(in1, 3)) max(P1(in1, 3))];
xbox2     = [min(P1(in1, 1)) max(P1(in1, 1))];
ybox2     = [min(P1(in1, 2)) max(P1(in1, 2))];
zbox2     = [min(P1(in1, 3)) max(P1(in1, 3))];
tri1      = (-delta+xbox1(1)<center1(:, 1)&center1(:, 1)<xbox1(2)+delta)&...
    (-delta+ybox1(1)<center1(:, 2)&center1(:, 2)<ybox1(2)+delta)&...
    (-delta+zbox1(1)<center1(:, 3)&center1(:, 3)<zbox1(2)+delta);
tri2      = (-delta+xbox2(1)<center2(:, 1)&center2(:, 1)<xbox2(2)+delta)&...
    (-delta+ybox2(1)<center2(:, 2)&center2(:, 2)<ybox2(2)+delta)&...
    (-delta+zbox2(1)<center2(:, 3)&center2(:, 3)<zbox2(2)+delta);
a = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
patch('Faces', t1(tri1, :), 'Vertices', P1, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c');
patch('Faces', t2(tri2, :), 'Vertices', P2, 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'FaceColor', 'y');
axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid off;
set(gcf,'color','w'); hold on;
plot3(P1(in1, 1), P1(in1, 2), P1(in1, 3), 'r.', 'MarkerSize', 15);
text(P1(in1, 1), P1(in1, 2), P1(in1, 3), num2str(in1), 'FontSize', 15);

title('Enclosed nodes: mesh #1 - cyan; mesh #2 - yellow');



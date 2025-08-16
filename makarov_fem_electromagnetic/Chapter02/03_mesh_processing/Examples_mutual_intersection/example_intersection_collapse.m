clear all; close all;
warning off; %  to avoid globals
%   SYNTAX 
%   example_intersection_collapse 
%   DESCRIPTION 
%   This script implements the mesh intersection algorithm for two meshes
%   and outputs difference mesh, union mesh, and intersection mesh
%   (performs basic Boolean operations). Additionally, this function
%   collapses small edges at the intersection boundary to improve mesh
%   quality. Both meshes A and B must be 2 manifolds
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%%  Mesh A (*.mat file)
FileName1 = 'VHP_Patella_left.mat';
load(FileName1);
P1 = P; t1 = t;
str.MasterMeshSize = size(t1, 1);
str.MasterMeshQ = min(simpqual(P1, t1));

%%  Mesh B (*.mat file)
FileName2 = 'VHP_Patella_left_shifted.mat';
load(FileName2);
P2 = P; t2 = t;
str.SlaveMeshSize  = size(t2, 1);
str.SlaveMeshQ = min(simpqual(P2, t2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Establish the necessary connectivity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
edges1      = meshconnee(t1);
edges_t1    = meshconnte(t1, edges1);
edges2      = meshconnee(t2);
edges_t2    = meshconnte(t2, edges2);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Establish a complete map of intersecting edges/triangles in steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%   Step 1 For every edge of mesh#2 find triangle(s) of mesh #1 intersected by
%   this edge and the corresponding intersection points (cell arrays si11, si12)
[si11, si12, si2] = meshedgeintersect(P1, t1, P2, edges2);
%   Step 2 For every edge of mesh#1 find triangle(s) of mesh #2 intersected by
%   this edge and the corresponding intersection points (cell arrays si11, si12)
[si21, si22, si1] = meshedgeintersect(P2, t2, P1, edges1);
%   Step 3 Collect all triangles of mesh #1 intersected by the edges of mesh #2
SI11 = unique(cell2mat(si11))';
%   Step 4 Collect all triangles of mesh #2 intersected by the edges of mesh #1
SI21 = unique(cell2mat(si21))';
%   Step 5 Collect all triangles of mesh #1 whose edges intersect triangles of mesh#2
AttachedTriangles = meshconnet(t1, edges1, 'manifold');
index = find(si1);
SI11extra = reshape(AttachedTriangles(index, :), 1, 2*length(index));
%   Step 6 Collect all triangles of mesh #2 whose edges intersect triangles of mesh#1
AttachedTriangles = meshconnet(t2, edges2, 'manifold');
index = find(si2);
SI21extra = reshape(AttachedTriangles(index, :), 1, 2*length(index));
%   Step 7 Collect all intersected or interesecting triangles of both meshes
SI11 = unique([SI11 SI11extra]);
SI21 = unique([SI21 SI21extra]);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Stop the script if iintersected or interesecting triangles were not found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(SI11)&isempty(SI21)
    display('The meshes have no intersections. The script will be stopped')
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Find extra segments/node pairs respecting all intersections in steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
st1  = size(t1, 1); st2  = size(t2, 1);
%  Step 1 For every triangle of mesh #1 find extra segments (node pairs) to be added
%  Duplicated points will be present
si13 = meshaddsegments(st1, st2, edges_t1, edges_t2, si11, si12, si21, si22, SI11, SI21);
%  Step 2 For every triangle of mesh #2 find extra segments (node pairs) to be added
%  Duplicated points will be present
si23 = meshaddsegments(st2, st1, edges_t2, edges_t1, si21, si22, si11, si12, SI21, SI11);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Create refined meshes without duplicated points in steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%  Create refined mesh #1 without duplicated points 
[tref1, Pref1] = meshrefined(t1, P1, si13);
str.RefinedMasterMeshSize = size(tref1, 1);
str.RefinedMasterMeshQ = min(simpqual(Pref1, tref1));
%  Create refined mesh #2 without duplicated points 
[tref2, Pref2] = meshrefined(t2, P2, si23);
str.RefinedSlaveMeshSize = size(tref2, 1);
str.RefinedSlaveMeshQ = min(simpqual(Pref2, tref2));
Pref10 = Pref1;
Pref20 = Pref2;
tref10 = tref1;
tref20 = tref2;
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Perform edge collapse at the boundary between two meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlength = 3; % This is a critical parameter
[Pref1, tref1, Pref2, tref2, steps] = meshcollapseboundary(Pref1, tref1, Pref2, tref2, minlength);
steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Perfom Boolean operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Center1 = (Pref1(tref1(:, 1),:) + Pref1(tref1(:, 2),:) + Pref1(tref1(:, 3),:))/3;         
Center2 = (Pref2(tref2(:, 1),:) + Pref2(tref2(:, 2),:) + Pref2(tref2(:, 3),:))/3;
[in1, trinumber1] = meshinoutstatus(Pref1, tref1, Center2);
[in2, trinumber2] = meshinoutstatus(Pref2, tref2, Center1);
%% A-B (difference of A and B)
Pm = [Pref1; Pref2];
tm = [tref1(~logical(in2), :); tref2(logical(in1), :)+size(Pref1, 1)];
[Pm, tm] = fixmesh(Pm, tm); 
%% A+B (union of A and B)
Pu = [Pref1; Pref2];
tu = [tref1(~logical(in2), :); tref2(~logical(in1), :)+size(Pref1, 1)];
[Pu, tu] = fixmesh(Pu, tu); 
%% AintB (intersection of A and B)
Pi = [Pref1; Pref2];
ti = [tref1(logical(in2), :); tref2(logical(in1), :)+size(Pref1, 1)];
[Pi, ti] = fixmesh(Pi, ti); 
%%  Mesh parameters
str.DifferenceSize      = size(tm, 1);
str.DifferenceQ         = min(simpqual(Pm, tm));
str.UnionSize           = size(tu, 1);
str.UnionQ              = min(simpqual(Pu, tu));
str.IntersectionSize    = size(ti, 1);
str.IntersectionQ       = min(simpqual(Pi, ti));
str
[nodes1, nodes2] = checknodes(Pref10, Pref20);
Nodes = Pref10(nodes1, :);
scrsz = get(0,'ScreenSize');

%%  Plot original meshes and all intersection nodes
a = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
patch('Faces', t1, 'Vertices', P1, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c');
patch('Faces', t2, 'Vertices', P2, 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'FaceColor', 'y');
axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid off;
set(gcf,'color','w'); hold on;
plot3(Nodes(:, 1), Nodes(:, 2), Nodes(:, 3), 'r.', 'MarkerSize', 25);
title('Intersection nodes');

%% Plot difference mesh and all intersection nodes
[nodes1, nodes2] = checknodes(Pref1, Pref2);
Nodes = Pref1(nodes1, :);
b = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
patch('Faces', tm, 'Vertices', Pm, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c'); 
xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid off;
set(gcf,'color','w'); hold on;
plot3(Nodes(:, 1), Nodes(:, 2), Nodes(:, 3), 'r.', 'MarkerSize', 25);
title('Difference mesh and intersection nodes');
axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
% P = Pm; t = tm;
% [P, t] = fixmesh(P, t);
% NewName =  strcat(FileName1(1:end-4), '_mod', '.mat');
% save(NewName, 'P', 't');

%% Plot union mesh and all intersection nodes
c = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
patch('Faces', tu, 'Vertices', Pu, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c');
axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid on;
set(gcf,'color','w'); hold on;
plot3(Nodes(:, 1), Nodes(:, 2), Nodes(:, 3), 'r.', 'MarkerSize', 25);
title('Union mesh and intersection nodes');
% P = Pu; t = tu;
% [P, t] = fixmesh(P, t);
% NewName =  strcat(FileName1(1:end-4), '_mod', '.mat');
% save(NewName, 'P', 't');

%% Plot intersection mesh and all intersection nodes
d = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid on;
set(gcf,'color','w'); hold on;
patch('Faces', ti, 'Vertices', Pi, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c'); hold on;
plot3(Nodes(:, 1), Nodes(:, 2), Nodes(:, 3), 'r.', 'MarkerSize', 25);
axis 'equal';    axis 'tight'; set(gca, 'YDir','normal');
title('Intersection mesh and intersection nodes');
% P = Pi; t = ti;
% [P, t] = fixmesh(P, t);
% NewName =  strcat(FileName1(1:end-4), '_mod', '.mat');
% save(NewName, 'P', 't');

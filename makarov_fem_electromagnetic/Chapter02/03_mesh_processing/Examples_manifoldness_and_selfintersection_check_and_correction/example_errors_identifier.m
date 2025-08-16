clear all; close all; warning off; 
%   SYNTAX 
%   example_errors_identifier 
%   DESCRIPTION 
%   This script identifies and removes all triangles attached to the
%   non-manifold edges and all self-intersecting triangles. The result is
%   a new mesh with holes, which is ready for healing (triangle addition or
%   further manual subtraction). The new mesh will be saved in ***_mod.mat
%   file. Multiple levels of triangle removal are possible.
%
%   Note: Given a meshes, which already has holes, the script will make
%   holes wider
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the mesh file for the self-intersection check');
load(FileName);
[str.MasterMeshQ str.MasterMeshIndex] = min(simpqual(P, t));
str

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Check the manifoldness condition first and collect suspicious triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges             = meshconnee(t);
AttachedTriangles = meshconnet(t, edges, 'nonmanifold');
NonManifoldAttached = [];
for m = 1:size(edges, 1)
    if length(AttachedTriangles{m})~=2        
        NonManifoldAttached = [NonManifoldAttached; AttachedTriangles{m}];
    end
end        
triangles = NonManifoldAttached;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Check self-intersections next and collect suspicious triangles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[si11, si12, si2] = meshedgeintersect(P, t, P, edges);
for m = 1:size(si11, 1)
    triangles = [triangles; si11{m}];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 1 Return if the mesh is good
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(triangles)
    display('The mesh is 2 manifold and has no self-intersections')
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 2 Remove all intersected/non-manifold triangles from the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t(triangles, :) = [];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 3 Remove all neighbor triangles - as many as necessary (control by M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
M = 0;
for m =1:M
    t(NonManifoldAttached, :) = [];
    [NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 4 Display the entire structure to see where the problem(s) occurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patch('Faces', t(:, :), 'Vertices', P, 'FaceColor', 'y', 'EdgeColor', 'k', 'FaceAlpha', 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 5 Remove unused nodes, display border triangles and border edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[P, t] = fixmesh(P, t);
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);
nodes = reshape(t(NonManifoldAttached, :), 1, 3*length(NonManifoldAttached));
nodes = unique(nodes);
for m = 1:size(edgesb, 1)
    marker1(m) = line('xdata', P(edgesb(m, :),1) ,'ydata', P(edgesb(m, :),2) ,'zdata', P(edgesb(m, :),3),...
    'color', 'b', 'linewidth', 4);           
end
patch('Faces', t(NonManifoldAttached, :), 'Vertices', P, 'FaceColor', 'c', 'EdgeColor', 'k', 'FaceAlpha', 1.0);
grid on; axis equal; axis tight; view(-136, -42)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Step 6 Save the mesh with holes for further processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
save(NewName, 'P', 't');

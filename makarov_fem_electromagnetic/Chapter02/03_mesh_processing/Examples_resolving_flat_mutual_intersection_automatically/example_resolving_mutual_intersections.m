clear all; close all;
warning off; %  to avoid globals
%   SYNTAX 
%   example_resolving mutual intersections
%   DESCRIPTION 
%   This prime script resolves mutual intersections for two meshes using local
%   topology-preserving mesh deformation in the direction of the unit
%   normal vector. The meshes must be 2-manifold. 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%%  Mesh A (*.mat file)
[FileName1, PathName] = uigetfile('*.mat','Select mesh #1 file for the intersection check');
load(FileName1);
P1 = P; t1 = t;
str.MasterMeshSize = size(t1, 1);

%%  Mesh B (*.mat file)
[FileName2, PathName] = uigetfile('*.mat','Select mesh #2 file for the intersection check');
load(FileName2);
P2 = P; t2 = t;
str.SlaveMeshSize  = size(t2, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Establish the necessary connectivity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
edges1      = meshconnee(t1);
%   Find triangles attached to every vertex (cell array)
s1 = meshconnvt(t1);
edges_t1    = meshconnte(t1, edges1);

edges2      = meshconnee(t2);
%   Find triangles attached to every vertex (cell array)
s2 = meshconnvt(t2);
edges_t2    = meshconnte(t2, edges2);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Main loop for moving the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 0;
alpha = -0.1;    % in mm
while 1 
    count = count + 1
    str.count = count;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Establish a complete map of intersecting edges/triangles in steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Stop the script if iintersected or interesecting triangles were not found
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(SI11)&isempty(SI21)
        display('The meshes have no intersections. The script will be stopped')               
        flag1       = checkintersection(P1, t1);
        flag2       = checkintersection(P2, t2);
        if ~(flag1+flag2)
            NewName =  strcat(FileName1(1:end-4), '_mod', '.mat');
            P = P1; t = t1;
            save(NewName, 'P', 't');
            NewName =  strcat(FileName2(1:end-4), '_mod', '.mat');
            P = P2; t = t2;
            save(NewName, 'P', 't');
        end
       %a = figure('Position',[1 scrsz(4)/4 scrsz(3)/1.5 scrsz(4)/1.5]);
       patch('Faces', t1, 'Vertices', P1, 'EdgeColor', 'b', 'FaceAlpha', 1.0, 'FaceColor', 'c');
       patch('Faces', t2, 'Vertices', P2, 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'FaceColor', 'y');
       axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
       xlabel('x, m'); ylabel('y, m'); zlabel('z, m'); view(-66, 19); grid off;        
       return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Move the nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %  Move nodes of mesh #1
    nodes1 = unique(t1(SI11, :));
    str.Nodes1 = length(nodes1);
    %   Find outer normal vectors
    normals = meshnormals(P1, t1);    
    %   Move selected nodes in the direction of the inner normal vector
    %   taking into account the local topology (important)
    index = nodes1;
    Ptemp = P1(index, :);
    for m = 1:length(index)
        n = index(m);
        averagenormal = sum(normals(s1{n}, :), 1)/length(s1{n});
        norm          = sqrt(dot(averagenormal, averagenormal));
        tangent       = norm;  
        Ptemp(m, :) = Ptemp(m, :) + alpha*(averagenormal/norm)/tangent^1.0;
    end 
    P1(index, :) = Ptemp;        
    
    %  Move nodes of mesh #2
    nodes2 = unique(t2(SI21, :));    
    str.Nodes2 = length(nodes2);
    %   Find outer normal vectors
    normals = meshnormals(P2, t2);    
    %   Move selected nodes in the direction of the inner normal vector
    %   taking into account the local topology (important)
    index = nodes2;
    Ptemp = P2(index, :);
    for m = 1:length(index)
        n = index(m);
        averagenormal = sum(normals(s2{n}, :), 1)/length(s2{n});
        norm          = sqrt(dot(averagenormal, averagenormal));
        tangent       = norm;  
        Ptemp(m, :) = Ptemp(m, :) + alpha*(averagenormal/norm)/tangent^1.0;
    end 
    P2(index, :) = Ptemp;       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str.MasterMeshQ     = min(simpqual(P1, t1));
    str.SlaveMeshQ      = min(simpqual(P2, t2));
    str
end



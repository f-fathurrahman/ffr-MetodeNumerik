clear all; close all; clc;
%   SYNTAX 
%   example_decimation_worst_triangle
%   DESCRIPTION 
%   This prime mesh-healing script sequentially decimates shortest edges of
%   the worst-quality triangles and tries to avoid self-intersections
%   and non-manifold meshes 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%   Load and check the original mesh
[FileName, PathName] = uigetfile('*.mat','Select the mesh file');
load(FileName);
[t, flag1]  = checkmanifold(t);
flag2       = checkintersection(P, t);
if flag1
    error('The mesh is not manifold');
end
if flag2
    error('The mesh has self-intersections');
end

%   Initial mesh data
str.tstart                  = size(t, 1);
str.Qstart                  = min(simpqual(P, t));
edges                       = meshconnee(t);
AttachedTriangles           = meshconnet(t, edges, 'manifold'); 
temp    = P(edges(:, 1),:) - P(edges(:, 2),:);
lengths = sqrt(dot(temp, temp, 2)); %   lengths of the edges    
[str.minlengthstart, dummy] = min(lengths);
str
    
%   Test edge lengths and perform edge collapse
M = 2; 
for m = 1:M
    m
    [str.MasterMeshQ str.MasterMeshIndex] = min(simpqual(P, t));
    %   Index into edges of the worst triangle
    index = meshconnte(t(str.MasterMeshIndex, :), edges);
    %   Perform edge collapse (try three different edge collapse schemes to
    %   avoid self-intersecting or non-manifold meshes)
    P0 = P; t0 = t;
    flag1 = zeros(3, 3);
    flag2 = zeros(3, 3);
    Qlocal= ones(3, 3)*str.MasterMeshQ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:3
        for j = 1:3
            [P, t] = meshcollapse(P0, t0, edges(index(i), :), AttachedTriangles(index(i), :), j);
            [t, flag1(i, j)] = checkmanifold(t);
            flag2(i, j)      = checkintersection(P, t);
            Qlocal(i, j)     = min(simpqual(P, t));
        end
    end
    %   Find the best case
    I = 0; J = 0; Q = str.MasterMeshQ;
    for i = 1:3
        for j = 1:3
            if (~flag1(i, j))&(~flag2(i, j))
                if Qlocal(i, j)>Q
                    I = i;
                    J = j;
                    Q = Qlocal(i, j)
                end
            end
        end
    end
    if (I*J)
        [P, t] = meshcollapse(P0, t0, edges(index(I), :), AttachedTriangles(index(I), :), J);
        save(strcat(FileName(1:end-4),'_mod'), 'P', 't');
    else
         error('Unable to proceed. Please fix the mesh');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edges                       = meshconnee(t);
    AttachedTriangles           = meshconnet(t, edges, 'manifold'); 
    temp                        = P(edges(:, 1),:) - P(edges(:, 2),:);
    lengths                     = sqrt(dot(temp, temp, 2)); %   the lengths of the edges    
    [str.minlength, index]      = min(lengths);    
    str.Qend                    = min(simpqual(P, t));
    str.tend                    = size(t, 1);
    str
end

str.Qend = min(simpqual(P, t));
str.tend = size(t, 1);
str










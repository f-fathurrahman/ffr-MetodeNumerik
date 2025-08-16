clear all; close all; clc;
%   SYNTAX 
%   example_decimation_shortest_edge
%   DESCRIPTION 
%   This script sequentially decimates the shortest edge of a mesh and attempts to avoid
%   self-intersections and non-manifold meshes
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
[str.minlengthstart, index] = min(lengths);
    
%   Test edge lengths and perform edge collapse
M = 5; %100;
for m = 1:M
    m
    %   Perform edge collapse (try three different edge collapse schemes to
    %   avoid self-intersecting or non-manifold meshes)
    P0 = P; t0 = t;
    for n = 1:3
        [P, t] = meshcollapse(P0, t0, edges(index, :), AttachedTriangles(index, :), n);
        [t, flag1] = checkmanifold(t);
        flag2      = checkintersection(P, t);
        if flag1
            display('The mesh is non-manifold');
        elseif flag2
            display('The mesh has self-intersections');
        else
            break;
        end
        if (flag1|flag2)&(n==3)
            error('Unable to proceed. Please fix the mesh');
        end
    end
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

save(strcat(FileName(1:end-4),'_mod'), 'P', 't');








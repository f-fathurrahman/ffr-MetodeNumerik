clear all; close all;
%   SYNTAX 
%   example_global_smoothing_lumped_laplace
%   DESCRIPTION 
%   This script performs global (or local) iterative Laplacian mesh
%   smoothing using the lumped-smoothing algorithm. The smoothing improves
%   triangle quality and increases min edge length but deforms the mesh
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the mesh file');
load(FileName);
[t, flag1]  = checkmanifold(t);
flag2       = checkintersection(P, t);
if flag1 error('The mesh is not manifold'); end
if flag2 error('The mesh has self-intersections'); end
%   Number of iteration steps (use large numbers)
M = 10;
%   Parameter alpha (use small numbers)
alpha = 0.1;

str.tri = size(t, 1);
str.Qstart = min(simpqual(P, t));
for m = 1:M
    m
    nodes = 1:size(P, 1);
    [P] = meshlaplace3Dlumped(P, t, nodes, alpha);
    str.Qend = min(simpqual(P, t));
end
str
[t, flag1]  = checkmanifold(t);
flag2       = checkintersection(P, t);
if flag1 error('The mesh is not manifold'); end
if flag2 error('The mesh has self-intersections'); end
save(strcat(FileName(1:end-4),'_mod'), 'P', 't');








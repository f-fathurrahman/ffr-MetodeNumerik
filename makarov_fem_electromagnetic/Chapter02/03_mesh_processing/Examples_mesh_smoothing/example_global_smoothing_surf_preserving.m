clear all; close all; 
%   SYNTAX 
%   example_global_smoothing_surf_preserving
%   DESCRIPTION 
%   This script performs global iterative surface-preserving Laplacian mesh
%   smoothing. The smoothing improves triangle quality and increases min
%   edge length but deforms the mesh
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

str.tri = size(t, 1);
str.Qstart = min(simpqual(P, t));
%   Number of iteration steps (use large numbers)
M = 1;
%   Parameters alpha, beta (use small numbers)
alpha = 0.1;
beta =  0.6;
[P, t] = meshlaplace3Dsp(P, t, alpha, beta, M);
str.Qend = min(simpqual(P, t));
str

[t, flag1]  = checkmanifold(t);
flag2       = checkintersection(P, t);
if flag1 error('The mesh is non-manifold'); end
if flag2 error('The mesh has self-intersections'); end
save(strcat(FileName(1:end-4),'_mod'), 'P', 't');








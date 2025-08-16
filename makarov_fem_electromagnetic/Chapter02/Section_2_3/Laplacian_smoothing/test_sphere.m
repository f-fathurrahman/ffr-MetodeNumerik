clear all;
%   SYNTAX
%   test_sphere
%   DESCRIPTION
%   This scripts evaluates function SPHERE
%   Note asymmetric deformation of the initially spherical domain!
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

R       = 1;         %   Sphere radius
Tr      = 500;       %   Approximate number of triangles
iter    = 20;        %   Number of iterations
[P, t] = sphere(R, Tr, iter);
size(t, 1)
clear all;
%   SYNTAX
%   test_circle
%   DESCRIPTION
%   This scripts evaluates function CIRCLE
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

R       = 1;         %   Circle radius
Tr      = 50;        %   Approximate number of triangles
iter    = 20;        %   Number of iterations
[P, t] = circle(R, Tr, iter);
size(t, 1)
clear all
%   SYNTAX
%   laplace3D
%   DESCRIPTION
%   This script performs 3D surface-preserving Laplacian mesh smoothing
%   based on existing Delaunay connectivity
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the object mesh file to open');
S = load(FileName, '-mat'); P = S.P; t = S.t; 

Quality         = min(simpqual(P, t))
graphics(P, t);  drawnow;


%   Algorithm
alpha   = 0.1; 
beta    = 0.6;
iter    = 5; 
for n = 1:iter
    [P, t] = meshlaplace3Dsp(P, t, alpha, beta, 1);
    clf; 
    graphics(P, t); title(strcat('Iter=', num2str(n)));
    Quality         = min(simpqual(P, t))
    pause(0.5);
end




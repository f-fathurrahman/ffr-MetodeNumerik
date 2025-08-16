clear all
%   SYNTAX
%   laplace3D
%   DESCRIPTION
%   This script performs 3D Laplacian mesh smoothing based on updated (through
%   MyRobustCrust)or existing Delaunay connectivity
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the object mesh file to open');
S = load(FileName, '-mat'); P = S.P; t = S.t; 

Quality         = min(simpqual(P, t))
graphics(P, t); drawnow;

n = 1;
while n<13
    step = n
    %   Triangles attached to every vertex (cell array)
    si = cell(size(P, 1), 1);
    for m = 1:size(P, 1)
        temp  = (t(:,1)==m)|(t(:,2)==m)|(t(:,3)==m);
        si{m} = find(temp>0);
    end        

    N = size(P, 1); 
    for m = 1:N                
        index = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
        index(index==m) = [];           %   Indexes into vertices connected to vertex m               
        Pnew(m, :) = sum(P(index, :))/length(index);   
    end

    % [t, normals]    = MyRobustCrust(Pnew);
    clf; graphics(Pnew, t); title(strcat('Iter=', num2str(n)));
    
    pause(0.5);
    Quality         = min(simpqual(Pnew, t))   
    P = Pnew;
    n = n + 1;
end
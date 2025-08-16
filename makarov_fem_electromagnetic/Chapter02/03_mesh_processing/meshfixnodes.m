function [P, t] = meshfixnodes(P, t)
%   SYNTAX
%   [P, t] = meshfixnodes(P, t)
%   DESCRIPTION
%   This function merges (nearly) coincident vertices, then removes
%   non-referenced vertices. It is a (slow) alternative to the corresponding
%   script of Per-Olof Persson. Tolerance parameter is given internally
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    tol = 1e-9;    
    %   Find coincident vertices
    if size(P, 2) == 3        
        sizex = max(abs(P(:, 1)));
        sizey = max(abs(P(:, 2)));
        sizez = max(abs(P(:, 3)));    
        dist = zeros(size(P, 1), size(P, 1));
        for m = 1:size(P, 1)
            distx = abs(P(:, 1) - P(m, 1))<tol*sizex;
            disty = abs(P(:, 2) - P(m, 2))<tol*sizey;
            distz = abs(P(:, 3) - P(m, 3))<tol*sizez;
            dist(m, :)  = distx&disty&distz;
        end
        dist = dist - diag(diag(dist));
    else
        %   Find coincident vertices
        sizex = max(abs(P(:, 1)));
        sizey = max(abs(P(:, 2)));
        dist = zeros(size(P, 1), size(P, 1));
        for m = 1:size(P, 1)
            distx = abs(P(:, 1) - P(m, 1))<tol*sizex;
            disty = abs(P(:, 2) - P(m, 2))<tol*sizey;
            dist(m, :)  = distx&disty;
        end
        dist = dist - diag(diag(dist));
    end            
    %   Replace coincident vertices in the array of triangles by a first
    %   nontrivial vertex
    %   Collect and sort non-referenced vertices
    Out = [];
    for m = 1:size(P, 1)
        index = find(dist(m, m:end)>0);      %   substitute vertex m for other coincident vertices
        index = index + m -1;
        Out = [Out index];
        for n = 1:length(index)
            t(t(:, 1)==index(n), 1) = m;
            t(t(:, 2)==index(n), 2) = m;
            t(t(:, 3)==index(n), 3) = m;
        end
    end
    Out = sort(unique(Out));
    %   Remove non-referenced vertices from the mesh entirely    
    Out(end + 1) = size(P, 1);
    for m = 1:length(Out)-1
        temp = t(:, 1)>Out(m)&t(:, 1)<Out(m+1);
        t(temp, 1) = t(temp, 1) - m;
        temp = t(:, 2)>Out(m)&t(:, 2)<Out(m+1);
        t(temp, 2) = t(temp, 2) - m;
        temp = t(:, 3)>Out(m)&t(:, 3)<Out(m+1);
        t(temp, 3) = t(temp, 3) - m;
    end
    P(Out(1: end-1), :) = [];
end
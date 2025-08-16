function [tMaster, PMaster] = meshrefined(t1, P1, si13)
%   SYNTAX
%   [t, P] = meshrefined(t1, P1, si13)
%   DESCRIPTION
%   This script is used in the intersection algorithm It creates refined
%   mesh #1 without duplicated points (adds extra triangles)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    tMaster = [];
    PMaster = [];
    for m = 1:size(t1, 1)
        if ~isempty(si13{m})    
            vector1 = P1(t1(m, 1), :) - P1(t1(m, 2), :); 
            vector2 = P1(t1(m, 1), :) - P1(t1(m, 3), :); 
            vector3 = cross(vector1, vector2);
            vector3 = cross(vector3, vector1); 
            xvector = vector1/norm(vector1);
            yvector = vector3/norm(vector3);
            basepoints = P1(t1(m, :), :); 
            [points, dummy, edges] = unique(si13{m}, 'rows', 'stable');  
            points = [points; basepoints];
            x = dot(points, repmat(xvector, size(points, 1), 1), 2);
            y = dot(points, repmat(yvector, size(points, 1), 1), 2);
            Edges = [];
            Edges(:, 1) = edges(1:2:end)';
            Edges(:, 2) = edges(2:2:end)';
            DT = delaunayTriangulation(x, y, Edges);
            Ptemp(:, 1) = x;
            Ptemp(:, 2) = y;
            Q = simpqual(Ptemp, DT.ConnectivityList);
            index = find(Q<2048*eps);
            tadd = DT.ConnectivityList; % local
            tadd(index, :) = [];
            Padd = points;              % local  
            clear Ptemp;                    
        else
            tadd = [1 2 3];             % local
            Padd = P1(t1(m, :), :);     % local              
        end
        tMaster = [tMaster; tadd+size(PMaster, 1)];
        PMaster = [PMaster; Padd];
    end
    [PMaster, tMaster] = fixmesh(PMaster, tMaster); 
end
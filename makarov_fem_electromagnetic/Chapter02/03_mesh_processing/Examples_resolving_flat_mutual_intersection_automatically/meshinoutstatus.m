function [in, trinumber] = meshinoutstatus(P, t, ObsPoints)
%   SYNTAX
%   [in, trinumber] = meshinoutstatus(P, t, ObsPoints)
%   DESCRIPTION
%   For an array ObsPoints, this function returns in = 1 if the point is
%   inside the manifold and zero otherwise
%   Inputs:
%   P - Vertex array of the mesh (N x 3)
%   t - Faces of the mesh (N x 3)
%   ObsPoints - array of points to be checked
%   Outputs:
%   The corresponding element of array in is 1 if the point is inside the
%   manifold and zero otherwise; 
%   trinumber - nearest triangle center for every point
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    % Find vertices of the triangles
    vert1 = P(t(:, 1),:);
    vert2 = P(t(:, 2),:);
    vert3 = P(t(:, 3),:);
    Center = (P(t(:, 1),:) + P(t(:, 2),:) + P(t(:, 3),:))/3;         % Calculating the center of the triangles
    in = zeros(size(ObsPoints, 1), 1);

    for m = 1:size(ObsPoints, 1)
        orig0 = ObsPoints(m, :);                                     % Finding the origin
        dist0 = Center - repmat(orig0, [size(t, 1), 1]);             
        dist  = dot(dist0, dist0, 2);                                % Finding distance to the nearest triangle
        [dummy, trinumber] = min(dist);
        dir0 = Center(trinumber, :) - orig0;
        dir0 = dir0/norm(dir0);                                      % Finding direction to the nearest triangle's center      
        T    = meshsegtrintersection(orig0, dir0, 1e6, P, t);
        T(T == 0) = [];                                              % Removing non-intersected triangles
        T = unique(T);                                               % Avoiding duplicate points
        if(rem(size(T, 1), 2)~=0)                                    % Checking for even or odd number of intersections
            in(m) = 1;                                               % If odd number of intersections, then the point is in
    end    
end

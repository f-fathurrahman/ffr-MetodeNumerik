function unitnormals = meshnormals(P, t)
%   SYNTAX
%   unitnormals = meshnormals(P, t)
%   DESCRIPTION
%   For a manifold mesh P, t, this function returns an array of outer
%   normal vectors for each triangle
%   Inputs:
%   P - Vertex array of the mesh (N x 3)
%   t - Facets of the mesh (N x 3)
%   Output:
%   unitnormals - Outer normal unit vectors for each triangle - N x 3
%   Authors: Vishal Rathi (vkrathi@wpi.edu)
%   Janakinadh Yanamadala (jyanamadala@wpi.edu), SNM (makarov@wpi.edu)
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

% Finding vertices of the triangle
vert1 = P(t(:, 1),:);
vert2 = P(t(:, 2),:);
vert3 = P(t(:, 3),:);
% Finding edges
edge1 = vert2 - vert1;
edge2 = vert3 - vert1;

normal = cross(edge1, edge2, 2);                                 % Calculating the normal of the triangle
length = sqrt(normal(:,1).^2+normal(:,2).^2+normal(:,3).^2);     % Calculating length of the normal
unitnormals= normal./(repmat(length,size(normal,3),3));          % Normalization of the normal vector
Center = (P(t(:, 1),:) + P(t(:, 2),:) + P(t(:, 3),:))/3;         % Calculating the center of the triangle

for m = 1:size(unitnormals, 1)
    orig0 = Center(m, :);                                        % Finding the origin
    dir0 = unitnormals(m, :);                                    % Finding the direction 
    T = meshsegtrintersection(orig0, dir0, 1e9, P, t);
    T(m) = 0;
    T(T ==0) = [];                                               % Avoiding origin points when t = 0
    T = unique(T);                                               % Avoiding duplicate points
    if(rem(size(T, 1), 2)~=0)                                    % Checking for even or odd number of intersections
        unitnormals(m,:)= -unitnormals(m,:);                     % If odd number of intersections, reverse the normal
    end
end
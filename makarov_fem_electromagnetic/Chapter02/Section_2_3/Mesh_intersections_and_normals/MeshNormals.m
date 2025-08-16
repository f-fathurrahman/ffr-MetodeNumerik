function unitnormals = MeshNormals(P, t)
%   SYNTAX
%   unitnormals = MeshNormals(P, t)
%   DESCRIPTION
%   For a manifold mesh P, t, this function returns an array of outer
%   normal vectors for each triangle
%   Inputs:
%   P - Vertex array of the mesh (N x 3)
%   t - Facets of the mesh (N x 3)
%   Output:
%   unitnormals - Outer normal unit vectors for each triangle - N x 3
%
%   For a more detailed description see a nearly equivalent function
%   meshnormals.m in the folder of major mesh functions
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

normal = cross(edge1, edge2,2);                                  % Calculating normal of the triangle
length = sqrt(normal(:,1).^2+normal(:,2).^2+normal(:,3).^2);     % Calculating length of the normal
unitnormals= normal./(repmat(length,size(normal,3),3));          % Normalization of the normal vector
Center = (P(t(:, 1),:) + P(t(:, 2),:) + P(t(:, 3),:))/3;         % Calculating the center of the triangle

for m = 1:size(unitnormals, 1)
    orig0 = Center(m, :);                                        % Finding origin
    dir0 = unitnormals(m,:);                                     % Finding direction
    orig  = repmat(orig0, size(vert1, 1),1);                     % Making dimensions similar
    dir  = repmat(dir0, size(vert1, 1),1);
    dist  = 1e10*ones(size(vert1, 1),1);
    [t] = SegmentTriangleIntersection(orig, dir, vert1, vert2, vert3,dist);
    t(m) = 0;
    t(t ==0) = [];                                               % Avoiding origin points when t = 0
    t = unique(t);                                               % Avoiding duplicate points
    if(rem(size(t,1),2)~=0)                                      % Checking for even or odd number of intersections
        unitnormals(m,:)= -unitnormals(m,:);                     % If odd number of intersections, reverse the normal
    end
end
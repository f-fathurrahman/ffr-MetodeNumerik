function d = SegmentTriangleIntersection(orig, dir, vert1, vert2, vert3, dist)
%   SYNTAX
%   d = SegmentTriangleIntersection(orig, dir, vert1, vert2, vert3, dist)
%   DESCRIPTION
%   This function checks whether or not a segment characterized by orig,
%   dir, dist intersects a triangle with vertices vert1, ver2, vert3
%   Input:
%   orig - Origin of the segment (N x 3)
%   dir - Normalized direction of the segment from origin (N x 3)
%   dist - Length of the segment (N x 1)
%   vert1, vert2, vert3 - Vertices of a single triangle or of all triangles
%   in a triangular mesh with N triangles
%   Output:
%   d -  Distances of points of intersection from the origin of the segment
%   (N x 1). If there is no intersection, the corresponding field is zero.
%   The tolerance is given internally
%
%   For a more detailed description see a nearly equivalent function
%   MESHSEGTRINTERSECTION.m in the folder of major mesh functions
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

% Initialization of u,v and d
u = zeros (size(vert1,1),1);
d = u; v = u;

% Finding edges
edge1 = vert2 - vert1;
edge2 = vert3 - vert1;

tvec = orig - vert1;            % Distance to vert1 from segment origin
pvec = cross(dir, edge2, 2);    % Parameter to calculate u
det  = dot(edge1, pvec, 2);     % Determinant of matrix M

parallel = abs(det)< eps;       % To test edges parallel with the segment
if all(parallel)                % If all parallel then no intersections
    return;
end

det(parallel) = 1;              % To avoid division by zero
inv_det = 1.0 ./ det;           % Find inverse of the determinant
u = dot(tvec,pvec,2);           % Calculate u parameter
u = u.*inv_det;

% Conditional tests for u and v
layer1 = (~ parallel & u<0 | u>1);
if all(layer1)
    return;
end

qvec (~layer1,:) = cross(tvec(~layer1,:), edge1(~layer1,:), 2);             % Parameter to calculate v
v (~layer1,:) = dot(dir(~layer1,:),qvec(~layer1,:),2).*inv_det(~layer1,:);  % Calculate v
layer2 = (v<=0 | u+v>1);
if all(layer2)
    return;
end

layer = (~layer1&~layer2);
d(layer,:) = dot(edge2(layer,:),qvec(layer,:),2).*inv_det(layer,:);         % Calculate d
d(d<0 | d>dist) = 0;                                                        % Compare distance and d
d(parallel) = 0;                                                            % Avoid any values of d in strictly parallel cases
d(isnan(d))= 0;                                                             % Avoiding NaN (Not-a-Number) when right-angled triangles are present
end



function [si11, si12, si2] = meshedgeintersect(P1, t1, P2, edges2);
%   SYNTAX 
%   [si11, si12, si2] = meshedgeintersect(P1, t1, P2, edges2) 
%   DESCRIPTION 
%   For every edge of mesh#2, this function returns triangle(s) of mesh #1
%   intersected by this edge and the corresponding intersection points
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    si11 = cell(size(edges2, 1), 1);    %   Intersected triangles
    si12 = cell(size(edges2, 1), 1);    %   Intersection points
    si2  = zeros(size(edges2, 1), 1);   %   Edge indicator for edges of mesh #2 
                                        %   intersecting mesh #1
%     %   Remove edges which are not inside of the convex hull of mesh#1
%     T       = delaunayTriangulation(P1);
%     e1      = pointLocation(T, P2(edges2(:, 1), :));
%     e2      = pointLocation(T, P2(edges2(:, 2), :));
%     index   = find(~(isnan(e1)&isnan(e2)));   %   a column
    index = [1:size(edges2, 1)]';
    %   Perform the loop for all other edges
    for n = 1:size(index, 1)
        m    = index(n);
        orig0= P2(edges2(m, 1), :);
        dest0= P2(edges2(m, 2), :);
        dir0 = dest0 - orig0;
        dist0= sqrt(dir0*dir0');
        dir0 = dir0/dist0;        
        t = meshsegtrintersection(orig0, dir0, (1-1e-9)*dist0, P1, t1);
        temp = find(abs(t)>1e-9*dist0);
        if ~ isempty(temp)
            si2(m) = 1;
            si11{m} = temp;    
            orig  = repmat(orig0, length(temp), 1);
            dir   = repmat(dir0, length(temp),  1);
            intersection = orig + repmat(t(temp), 1, 3).*dir;    %  Finding intersection points
            si12{m} = intersection; 
        end
    end
end
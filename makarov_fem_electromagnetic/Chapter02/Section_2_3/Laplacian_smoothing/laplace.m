function [Pnew, t] = laplace(P, I, type)
%   SYNTAX 
%   [Pnew, t] = laplace(P, I, type)
%   DESCRIPTION 
%   This function implements 2D and 3D iterative Laplacian smoothing with
%   retriangulation and fixed boundary (only in 2D)
%   Inputs:
%   P - array of vertices
%   I - first I vertices will not be moved - fixed vertices
%   type = 1 for standard Laplacian smoothing
%   type = 2 for lumped Laplacian smoothing
%   type = 3 for Centroid Voronoi Tessellation (CVT) smoothing  (t-centers)
%   type = 4 for Weighted Centroid of Circumcenters (WCC) smoothing  (c-centers)
%   Outputs: 
%   Pnew - array of nodes
%   t - triangulation
%  
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    %   Unconstrained Delaunay triangulation (2D or 3D)
    DT  = delaunayTriangulation(P);
    t   =  DT.ConnectivityList;
    
     if size(P, 2) == 2 % 2D
        DT  = delaunayTriangulation(P);
        t   =  DT.ConnectivityList;
     else  % 3D
        DT  = delaunayTriangulation(P);       
        t   = freeBoundary(DT);
        DT   = triangulation(t, P);                 
     end
    
    si  = vertexAttachments(DT);    %   Triangles attached to every vertex (cell array)
    N   = size(P, 1);
    
    %   Areas   
    d12 = P(t(:,2),:)-P(t(:,1),:);
    d13 = P(t(:,3),:)-P(t(:,1),:);
    A   = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;

    Pnew = P;
    switch type
        case 1  %   Laplacian smoothing      
            for m = I+1:N                
                index = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
                index(index==m) = [];           %   Indexes into vertices connected to vertex m               
                Pnew(m, :) = sum(P(index, :))/length(index);   
            end
        case 2  %   lumped Laplacian smoothing
            for m = I+1:N
                index = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
                index(index==m) = [];           %   Indexes into vertices connected to vertex m
                Pnew(m, :) = 2/3*sum(P(index, :))/length(index) + 1/3*P(m, :);
            end       
        case 3  %   Centroid Voronoi Tessellation (CVT) smoothing
            ic      = incenter(DT);            %   Centroids of attached triangles
            for m = I+1:N
                tcenters = ic(si{m}', :);                             
                tareas   = A(si{m}');
                Pnew(m, :) = sum(tcenters.*repmat(tareas, 1, 2), 1)/sum(tareas);
            end       
        case 4  %    Weighted Centroid of Circumcenters (WCC) smoothing
            cc      = circumcenter(DT);        %   Circumcenters of attached triangles
            for m = I+1:N
                ccenters = cc(si{m}', :);                  
                tareas   = A(si{m}');           
                Pnew(m, :) = sum(ccenters.*repmat(tareas, 1, 2), 1)/sum(tareas);
            end      
    end
%     %   Shrinkage avoidance (return the vertices back to the boundary for the sphere)
%     if R > 0
%         for m = I+1:N
%             Pnew(m, :) = R*Pnew(m, :)/sqrt(sum(Pnew(m, :).*Pnew(m, :)));
%         end
%     end    
end
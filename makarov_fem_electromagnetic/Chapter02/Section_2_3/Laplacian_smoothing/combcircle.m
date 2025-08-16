function [P, t] = combcircle(R, Tri, iter, strint)
%   SYNTAX
%   [P, t] = combcircle(R, Tri, iter, strint)
%   DESCRIPTION
%   This function creates a mesh for a base circle of radius R, with
%   approximately Tri triangles, and with circular subdomains desribed in
%   structure strint using Laplacian smoothing with retriangulation
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    %%   Main circle
    %   Determine the required number of subdivisions along circumference
    N = ceil(sqrt(Tri)/0.4);
    
    %   Create a lattice for equilateral triangles in the enclosing rectangle
    sizex   = 2*pi*R/N;             
    sizey   = sizex*tan(pi/3)/2; 
    x       = -R-sizex/2:sizex:R+sizex/2; 
    y       = -R-sizey/2:sizey:R+sizey/2;
    
    [X, Y]  = meshgrid(x, y);
    X(1:2:end, :) = X(1:2:end, :) + sizex/4;     %   Shift odd rows
    X(2:2:end, :) = X(2:2:end, :) - sizex/4;     %   Shift even rows
    P = [X(:), Y(:)];                            %   List of vertices 
    
    %   Remove vertices outside the circle and close to its boundary
    temp        = dot(P, P, 2)>(R-0.5*sizex)^2;
    P(temp, :)  = [];

    %   Add boundary nodes and boundary connectivity for the main (outer) boundary
    [Pb, vert] = combboundary(N, R, 0, 0, [], []);
    
    %   Add boundary nodes and boundary connectivity for internal boundaries
    for m = 1:size(strint.R, 2)
       [Pb, vert] = combboundary(strint.N(m), strint.R(m), strint.x(m), strint.y(m), Pb, vert);
    end
 
    %%  Combine meshes
    P    = [Pb; P];
    
    %   Create and view the initial structured mesh (constrained Delaunay)
    DT  = delaunayTriangulation(P, vert);
    t   =  DT.ConnectivityList;

    %   Graphics
    fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
    t = combdomain(P, t, strint);
    for n = 1:size(strint.R, 2)
        fv.faces = t(t(:, 4)==n, 1:3); fv.vertices = P; patch(fv, 'FaceColor', strint.c(n));            
    end     
    axis equal; view(0, 90); grid on; str.iter = 0; str.minquality = min(simpqual(P, t)); str    
    disp('Press Enter'); 
    pause
    
    %   Perform Laplacian smoothing iteratively
    for m = 1:iter        
        [dummy, index] = min(simpqual(P, t));  
        %   Remove the lowest-quality triangle (delete its node if not on the boundary)     
        if m<iter
            node = t(index, find(t(index, :)>size(Pb, 1)));
            if length(node)==1 P(node, :) = []; end;
        end
        [P, t] = laplace(P, size(Pb, 1), 4);     
        DT  = delaunayTriangulation(P, vert);
        t   =  DT.ConnectivityList;
        %   Graphics          
        clf; fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
        t = combdomain(P, t, strint);
        for n = 1:size(strint.R, 2)
            fv.faces = t(t(:, 4)==n, 1:3); fv.vertices = P; patch(fv, 'FaceColor', strint.c(n));            
        end          
        t(:, 4) = 0;
        axis equal; view(0, 90); grid on; xlabel('x'); ylabel('y'); drawnow;
        str.iter = m; str.minquality = dummy; str
        disp('Press Enter');
        pause
    end
    
end

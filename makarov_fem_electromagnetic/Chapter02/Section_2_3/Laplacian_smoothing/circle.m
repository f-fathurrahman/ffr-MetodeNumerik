function [P, t] = circle(R, Tri, iter)
%   SYNTAX
%   [P, t] = circle(R, Tri, iter)
%   DESCRIPTION
%   This function creates a mesh for a base circle of the radius R and
%   with approximately Tri triangles using Laplacian smoothing
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

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

    %   Add boundary nodes    
    phi = [0:2*pi/N:2*pi*(N-1)/N];
    bx = R*cos(phi); by = R*sin(phi);
    I = size(bx, 2);
    P = [[bx; by]'; P]; 
 
    %   Create and view the initial structured mesh 
    DT  = delaunayTriangulation(P);
    t   =  DT.ConnectivityList;
    
       
    fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
    axis equal; view(0, 90); grid on;
    str.iter = 0; str.quality = min(simpqual(P, t)); str
    disp('Press Enter'); 
    pause
  
    %   Perform Laplacian smoothing
    for m = 1:iter
        [P, t] = laplace(P, I, 4);        
        clf; fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y');  
        axis equal; view(0, 90); grid on; xlabel('x'); ylabel('y'); drawnow;
        str.iter = m; str.quality = min(simpqual(P, t)); str
        disp('Press Enter'); 
        pause
    end
 
end

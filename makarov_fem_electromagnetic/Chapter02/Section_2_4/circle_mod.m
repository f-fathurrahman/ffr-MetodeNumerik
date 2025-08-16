function [P, t, I] = circle_mod(R, Tri)
%   SYNTAX
%   [P, t, I] = circle_mod(R, Tri)
%   DESCRIPTION
%   This function creates a mesh for a base circle of the radius R and with
%   approximately Tri triangles using Laplacian smoothing. I is the number
%   of boundary nodes placed first
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
       
    %   Perform Laplacian smoothing with up to 20 iterations
    iter = 20;
    for m = 1:iter
        q1 = min(simpqual(P, t));
        [P, t] = laplace_mod(P, I, 4);  
        if  min(simpqual(P, t))<q1
            break;
        end
    end
end

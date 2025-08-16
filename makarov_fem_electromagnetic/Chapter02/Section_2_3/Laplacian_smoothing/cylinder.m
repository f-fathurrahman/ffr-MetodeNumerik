function [P, t] = cylinder(R, H, Tri);
%   SYNTAX
%   [P, t] = cylinder(R, H, Tri);
%   DESCRIPTION
%   This function creates a mesh for a base cylinder with radius R, height
%   H, and orientation along the z-axis. Tri is approximate number of
%   triangles using Laplacian smoothing for cylinder caps
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
    
    %   Generate cylinder caps first
    %   Determine the required number of subdivisions along circumference
    N = ceil(sqrt(Tri/(0.32 + H/(pi*R))));
   
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
 
    %   Perform Laplacian smoothing
    iter = 10;
    for m = 1:iter
        [P, t] = laplace(P, I, 4);        
    end
    
    %  Determine z-divisions
    blength = sqrt(R^2 + R^2 -2*R^2*cos(2*pi/N));       %   boundary edge length
    hsize   = blength;                                  %   anticipated size of the z-divisions
    NN      = ceil(H/hsize);                            %   number of divisions along the height
    hsize   = H/NN;                                     %   true size of the z-divisions
    
    %   Replicate the circle mesh at all different z-locations
    P(:, 3) = - H/2; Pcap = P; pcap = size(Pcap, 1);
    for m = 1:NN
        P(end+1:end+pcap, :) = Pcap;
        P(end+1-pcap:end, 3) = -H/2 + hsize*m;
    end
    
    %   Remove all unnecessary inner vertices
    tol = 0.01;
    temp        = (dot(P(:, 1:2), P(:, 1:2), 2)<(R-tol*blength)^2)&abs(P(:, 3))<H/2-tol*hsize;
    P(temp, :)  = []; 
  
    %   Unconstrained Delaunay triangulation (2D and 3D)
    DT  = delaunayTriangulation(P);       
    t   = freeBoundary(DT);
end

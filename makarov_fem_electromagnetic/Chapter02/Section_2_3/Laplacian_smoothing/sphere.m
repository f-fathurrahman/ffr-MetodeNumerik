function [P, t] = sphere(R, Tri, iter)
%   SYNTAX
%   [P, t] = sphere(R, Tri, iter)
%   DESCRIPTION
%   This function creates a mesh for a sphere of radius R centered at
%   origin using iterative Laplacian smoothing. The initial nodal points
%   are assigned exactly on the sphere surface. Tri is approximate number
%   of triangles
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
 
    %   Create the initial structured grid - a series of circles centered at the z-axis
    zlength = 2.5*R/sqrt(Tri);          %   Circle spacing at equator
    Nz      = 2*R/zlength;              %   Divisions along  the z-axis
    z       = -R*cos(pi*[0:Nz]/Nz);     %   Non-uniform grid along the z-axis 
                                        %   (simplified, but keeps two asymptotics)
    r = sqrt(R^2 - z.^2);               %   Circle radii    
    temp = zlength/(tan(pi/3)/2)*2;
    Ptemp = [];                         %   Accumulate the array of vertices
    for m = 1:Nz
        N = floor(2*pi*r(m)/temp);
        phi = [0:2*pi/N:2*pi*(N-1)/N] + mod(m, 2)*pi/N;    %   interleaving
        bx = r(m)*cos(phi); by = r(m)*sin(phi);
        Ptemp = [Ptemp  [bx; by; z(m)*ones(size(bx))] ];
    end
    P   = Ptemp'; I = size(P, 1);
    
    %   Create and view the initial structured mesh 
    DT  = delaunayTriangulation(P);       
    t   = freeBoundary(DT);     % Keep only the surface mesh
       
    fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
    axis equal; view(160, 60); grid on;
    str.iter = 0; str.quality = min(simpqual(P, t)); str
    disp('Press Enter'); 
    pause
    
    %   Perform Laplacian smoothing
    for m = 1:iter
        [P, t] = laplace(P, 0, 1);        
        clf; fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y');  
        axis equal; view(160, 60); grid on; drawnow;
        xlabel('x'); ylabel('y');
        str.iter = m; str.quality = min(simpqual(P, t)); str
        disp('Press Enter'); 
        pause
    end
end
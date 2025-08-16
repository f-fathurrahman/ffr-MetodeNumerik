function [P, t, I] = plate_mod(L, W, Tri, iter);
%   SYNTAX
%   [P, t, I] = plate_mod(L, W, Tri, iter);
%   DESCRIPTION
%   This function creates a mesh for a base plate of length L and width W
%   in the xy plane using iterative Laplacian smoothing. The initial nodal
%   points are assigned in a triangular lattice. Tri is approximate number
%   of triangles. I is the number of boundary nodes placed first
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
    
    %   Create a lattice for equilateral triangles within the enclosing rectangle
    N = ceil(sqrt(Tri*L/(2.2*W))); sizex = L/(N-1); sizey  = sizex*sqrt(3)/2; 
    x       = -L/2+sizex/2:sizex:L/2-sizex/2; 
    y       = -W/2+sizey/2:sizey:W/2-sizey/2;
    [X, Y]  = meshgrid(x, y);
    X(1:2:end, :) = X(1:2:end, :) + sizex/4;     %   Shift odd rows
    X(2:2:end, :) = X(2:2:end, :) - sizex/4;     %   Shift even rows
    P = [X(:), Y(:)];                            %   List of inner vertices 
    
    %   Add boundary nodes (clockwise)
    x = linspace(-L/2, L/2, ceil(L/sizex));
    y = linspace(-W/2, W/2, ceil(W/sizey));
    temp1 = [-L/2*ones(1, length(y));   y                ];
    temp2 = [x(2:end-1);    +W/2*ones(1, length(x)-2)    ];
    temp3 = [+L/2*ones(1, length(y));   y(end:-1:1)      ];
    temp4 = [x(end-1:-1:2);     -W/2*ones(1, length(x)-2)]; 
    temp  = [temp1 temp2 temp3 temp4]';
    I = size(temp, 1);
    P = [temp; P]; 
    
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
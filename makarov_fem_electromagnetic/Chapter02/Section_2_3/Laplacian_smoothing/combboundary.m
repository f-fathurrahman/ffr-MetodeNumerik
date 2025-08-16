function [Pnew, vertnew] = combboundary(N, R, x, y, P, vert)
%   SYNTAX
%   [Pnew, vertnew] = combboundary(N, R, x, y, P, vert)
%   DESCRIPTION
%   This function adds boundary nodes and boundary connectivity for
%   internal circles
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    % Add boundary nodes
    phi     = [0:2*pi/N:2*pi*(N-1)/N];
    bx      = R*cos(phi) - x; 
    by      = R*sin(phi) - y;
    I0      = size(bx, 2);
    P0      = [bx; by]'; 
    Pnew    = [P; P0];
    
    %   Add boundary vertixes (connectivity)
    vert0       = [[1:I0]; [2:I0 1]]';
    vertnew     = [vert; vert0+size(vert, 1)];

end


function [P, t] = brick0(L, W, H, Nx, Ny, Nz, indicator);
%   SYNTAX
%   [P, t] = [P, t] = brick0(L, W, H, Nx, Ny, Nz, indicator);
%   DESCRIPTION
%   This function creates a structured mesh for a base brick on
%   the size L(length-x) by W(width-y) by H(height-z). Use indicator=1 for
%   a non-uniform mesh.
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    %   Create six plates 
    [P1, t1] = plate0(L, W, Nx, Ny, indicator, 0, 0, 1);
    P1       = move(P1, 0, 0, -H/2);                  % bottom
    [P2, t2] = plate0(L, W, Nx, Ny, indicator, 0, 0, 1);
    P2       = move(P2, 0, 0, +H/2);                  % top
    [P3, t3] = plate0(W, H, Ny, Nz, indicator, 0, 0, 1);
    P3       = rotatex(P3, 90);    
    P3       = rotatez(P3, 90);    
    P3       = move(P3, -L/2, 0, 0);                  % left
    [P4, t4] = plate0(W, H, Ny, Nz, indicator, 0, 0, 1);
    P4       = rotatex(P4, 90);    
    P4       = rotatez(P4, 90);     
    P4       = move(P4, +L/2, 0, 0);                  % right      
    [P5, t5] = plate0(L, H, Nx, Nz, indicator, 0, 0, 1);
    P5       = rotatex(P5, 90);    
    P5       = move(P5, 0, -W/2, 0);                  % rear
    [P6, t6] = plate0(L, H, Nx, Nz, indicator, 0, 0, 1);
    P6       = rotatex(P6, 90);    
    P6       = move(P6, 0, +W/2, 0);                  % front 
    %   Unite mesh
    [P, t] = meshcombine(P1, P2, P3, P4, P5, P6, t1, t2, t3, t4, t5, t6);
    %   Merge coincident vertrices (if any)
    [P, t] = meshfixnodes(P, t);
end

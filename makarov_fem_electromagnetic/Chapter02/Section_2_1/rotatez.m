function Pnew = rotatez(P, zangle)
%   SYNTAX
%   P = rotatez(P, zangle) where the angle is in degrees
%   DESCRIPTION
%   This function rotates the mesh about the z-axis
%   The local coordinate system may be used with the origin at the center of gravity
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    anglez = zangle/180*pi;
    %   Rotation about the z-axis
    LCX  = mean(P(:, 1)); LCX = 0;
    LCY  = mean(P(:, 2)); LCY = 0;      
    Pnew(:, 1)  = +(P(:, 1) - LCX)*cos(anglez) - (P(:, 2) - LCY)*sin(anglez);
    Pnew(:, 2)  = +(P(:, 1) - LCX)*sin(anglez) + (P(:, 2) - LCY)*cos(anglez);
    Pnew(:, 3)  = +P(:, 3); 
    Pnew(:, 1) = Pnew(:, 1) + LCX;
    Pnew(:, 2) = Pnew(:, 2) + LCY;
end


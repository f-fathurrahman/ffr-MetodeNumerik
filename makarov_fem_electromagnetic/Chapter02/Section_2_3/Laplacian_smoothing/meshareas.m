function A = meshareas(P, t) 
%   SYNTAX 
%   A = meshareas(P, t);
%   DESCRIPTION 
%   This function returns areas of all triangles in the mesh in array A
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    d12     = P(t(:,2),:)-P(t(:,1),:);
    d13     = P(t(:,3),:)-P(t(:,1),:);
    temp    = cross(d12, d13, 2);
    A       = 0.5*sqrt(dot(temp, temp, 2));
end
   
    

clear all;
%   SYNTAX
%   refine2_plate
%   DESCRIPTION
%   This script models adaptive mesh refinement for a plate - performs
%   subdivision of triangles with the larger error. Laplacian smoothing of
%   the resulting mesh is used 
%   L, W - plate length and width
%   Tri - approximate number of triangles in the initial nearly uniform mesh 
%   P - input array of nodes
%   t - input connectivity
%   I - number of boundary nodes - placed first in array P
%   error - error value for each triangle
%   percentage - triangle percentage to be refined
%   iter - number of mesh refinement steps
%   par:
%   par = 0->Only edges adjacent to two triangles in question are refined
%   par = 1->All edges of the triangle in question are refined
%   method/type - Laplacian method (1-4)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%   Input data
L           = 1; 
W           = 1;
Tri         = 200;
percentage  = 5;
iter        = 12;
par         = 1;
method      = 4;

%   Initial mesh
[P, t, I] = plate_mod(L, W, Tri);
fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y');     
axis equal; view(0, 90); grid on; str.iter = 0; str.minquality = min(simpqual(P, t)); str.triangles = size(t, 1); str    
disp('Press ENTER'); pause

for n = 1:iter    
    %   triangle centers and areas
    TriCenter = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) +  P(t(:, 3), :));
    d12 = P(t(:,2),:)-P(t(:,1),:);
    d13 = P(t(:,3),:)-P(t(:,1),:);
    A   = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;
    %   definition of the local error
    d = drectangle(TriCenter, L, W); d(d>0) = 0; 
    error = 1./sqrt(abs(d)).*A.^0.5;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf;
    subplot(1, 2, 1);
    p = patch('faces', t, 'vertices', P); cdata =error; set(p,'FaceColor','flat','FaceVertexCData',cdata, 'CDataMapping','scaled');
    axis equal;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [P, t, I] = plate_refine(P, t, I, error, percentage, par, method); %   Adding nodes and refinement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1, 2, 2); patch('faces', t, 'vertices', P, 'FaceColor', 'y');
    axis equal;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str.iter = n; str.minquality = min(simpqual(P, t)); str.triangles = size(t, 1); str     
    drawnow; disp('Press ENTER'); pause;
end

dt                  = triangulation(t(:, 1:3), P);
[IC, RIC]           = incenter(dt);
trianglesizeratio   = max(RIC)/min(RIC)


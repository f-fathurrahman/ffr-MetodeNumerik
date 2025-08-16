clear all;
%   Create a uniform mesh for an ellipsoid with approximately Tr triangles
%   Original code: DISTMESH 2004-2012 Per-Olof Persson
%   Adopted by SNM Summer 2012 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

R = 0.1; Tr = 300;

%   Distance function
fd = inline('sqrt(sum((P./repmat([1 1.5 1],length(P),1)).^2, 2))-R', 'P', 'R');

h0    = 5.5*R/sqrt(Tr);                             %   Approximate desired edge length
xmin  = -R; xmax = +R; ymin = -R; ymax = +R;        %   3D Bounding box
zmin  = -R; zmax = +R; 

dptol   = 0.005;                                    %   Termination criterion
ttol    = 0.1;                                      %   Criterion for retriangulation
Fscale  = 1.2;                                      %   Inner "internal pressure"
deltat  = 0.2;                                      %   Ratio between force and movement
geps    = 0.001*h0;                                 %   Relative boundary tolerance
deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation

%   Create initial uniform distribution in bounding box (MATLAB isosurface from grid)
[x, y, z] = ndgrid(xmin:h0:xmax, ymin:h0:ymax, zmin:h0:zmax);
pv        = isosurface(x, y, z, reshape(fd([x(:), y(:), z(:)], R), size(x)), 0);
P         = pv.vertices;

N    = size(P, 1);                                  %  Number of vertices N
Pold = inf;                                         %  For first iteration

count = 0;
while 1
    %   Retriangulation by the Delaunay algorithm
    if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
        Pold    = P;                                %   Save current positions
        DT  = DelaunayTri(P);                       %   3D Delaunay triangulation
        t   = freeBoundary(DT);                     %   3D Triangular mesh
        bars    = edges(TriRep(t, P));              %   All sorted mesh edges [M, 2]
        M       = size(bars, 1);                    %   Length of the array of edges
        count = count + 1;  
    end
    
    %   Move mesh points based on edge lengths L1 and forces F
    barvec    = P(bars(:, 1),:)-P(bars(:, 2),:);            %   List of edge vectors [M, 2]
    L1         = sqrt(sum(barvec.^2,2));                    %   Edge lengths [M, 1]
    L0        = Fscale*sqrt(sum(L1.^2)/M)*ones(M, 1);       %   Desired (average) edge lengths [M, 1]
    F         = max(L0-L1, 0);                              %   Edge forces [M, 1]
    Fvec      = repmat(F./L1, 1, 3).*barvec;                %   Edge forces (x,y,z components) [M, 3]
    Ftot      = zeros(N, 3);                                %   Total force for the vertices [N, 3]
    
    %   Add up force contributions for vertices (the original code - one line - runs much faster)
    %     for m = 1:N                                             
    %         ind1 = find(bars(:, 1)==m);
    %         ind2 = find(bars(:, 2)==m);
    %         Ftot(m, :) = sum(Fvec(ind1, :), 1) - sum(Fvec(ind2, :), 1);
    %     end    
    Ftot = full(sparse(bars(:,[1,1,1,2,2,2]), ones(size(F))*[1,2,3,1,2,3], [Fvec,-Fvec], N, 3));
    P = P + deltat*Ftot;                            %   Update node positions
    
    %   Bring all vertices back to the boundary
    dist  = feval(fd, P, R);                        %   Find distances
    dgradx = (feval(fd, [P(:, 1)+deps, P(:, 2), P(:, 3)], R)-dist)/deps;   %   Numerical
    dgrady = (feval(fd, [P(:, 1), P(:, 2)+deps, P(:, 3)], R)-dist)/deps;   %   gradient
    dgradz = (feval(fd, [P(:, 1), P(:, 2), P(:, 3)+deps], R)-dist)/deps;   %   gradient
    dgrad2 = dgradx.^2 + dgrady.^2 + dgradz.^2;
    P = P - [dist.*dgradx./dgrad2, dist.*dgrady./dgrad2, dist.*dgradz./dgrad2];     
                                                                           %   Projection
    %   Termination criterion: All interior nodes move less than dptol (scaled)
    factor = max(sqrt(sum((P-Pold).^2, 2))/h0);
    if factor<dptol, break; end
    if count>2^10 break; end;
end

count
cla; patch('vertices', P, 'faces',t,'edgecol','k','facecol',[.8, .9, 1]); view(0, 90); axis equal;
xlabel('x'); ylabel('y');
Triangles = size(t, 1)
quality = min(simpqual(P, t))
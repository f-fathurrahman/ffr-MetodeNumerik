clear all;
%   Create a uniform or nonuniform mesh  for a base plate (rectangle) with 
%   length L and width W, and with approximately Tr triangles
%   Original code: DISTMESH 2004-2012 Per-Olof Persson
%   Adopted by SNM Summer 2012 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

Tr  = 1000;             %   Approximate number of triangles
L   = 1;                %   Rectangle length
W   = 1;                %   Rectangle width
par = 0.2;              %   Mesh uniformity parameter (0-1) 

h    = waitbar(0.5, 'Please wait - computing triangular mesh');

h0    = 1.40*sqrt(L*W)/sqrt(Tr);                    %   Approximate desired edge length
xmin  = -L/2; xmax = +L/2; ymin = -W/2; ymax = +W/2;%   Bounding box

dptol   = 0.001;                                    %   Termination criterion
ttol    = 0.1;                                      %   Criterion for retriangulation  
Fscale  = 1.2;                                      %   Inner "internal pressure" (original)
if Tr<=200 & par<0.3
    Fscale = 1.10;
end
deltat  = 0.2;                                      %   Ratio between force and movement
geps    = 0.001*h0;                                 %   Relative boundary tolerance (original)
if par ==1
    geps    = 0.1*h0;                               %   Relative boundary tolerance (relaxed here)
end  
deps    = sqrt(eps)*h0;                             %   Accuracy of gradient calculation

%   Create initial distribution in bounding box 
%   Equilateral triangles - faster convergence
[x, y] = meshgrid(xmin:h0:xmax, ymin:h0*sqrt(3)/2:ymax);
x(2:2:end, :) = x(2:2:end, :) + h0/2;                 %   Shift even rows
P = [x(:), y(:)];                                     %   List of node coordinates

%   Remove vertices outside the region
dist = drectangle(P, L, W);                         %  A column of signed distances
P    = P(dist<geps, :);                             %  Keep only d<0 (geps) points  
N    = size(P, 1);                                  %  Number of vertices N
Pold = inf;                                         %  For first iteration

count   = 0;    
while 1      
    %   Retriangulation by the Delaunay algorithm
    if max(sqrt(sum((P-Pold).^2,2))/h0)>ttol        %   Any large movement?
        Pold    = P;                                %   Save current positions
        DT  = delaunayTriangulation(P);             %   2D Delaunay triangulation
        t   =  DT.ConnectivityList;                 %   List of triangles
        ic      = incenter(DT);                     %   Centroids of triangles
        dist    = drectangle(ic, L, W);             %   A column of signed distances
        t       = t(dist<-geps, :);                 %   Keep interior trianges
        bars    = edges(triangulation(t, P));       %   All sorted mesh edges [M, 2]
        M       = size(bars, 1);                    %   Length of the array of edges        
    end
    count = count + 1;     

    %   Move mesh points based on edge lengths L and forces F
    barvec  = P(bars(:, 1),:)-P(bars(:, 2),:);              %   List of edge vectors [M, 2]
    L1      = sqrt(sum(barvec.^2,2));                       %   Edge lengths [M, 1]       
    barsc   = 0.5*(P(bars(:, 1), :) + P(bars(:, 2), :));    %   Edge centers [M, 2]
    hbars   = hrectangle(barsc, L, W, par);                 %   Element sizes for edge centers [M, 1]
    hbars   = hbars/sqrt(sum(hbars.^2)/M);                  %   Element sizes normalized by its average [M, 1] 
    L0      = Fscale*sqrt(sum(L1.^2)/M)*hbars;              %   Desired (non-uniform) edge lengths [M, 1]        
    F       = max(L0-L1, 0);                                %   Edge forces [M, 1]
    Fvec    = repmat(F./L1, 1, 2).*barvec;                  %   Edge forces (x,y components) [M, 2]

    Ftot = full(sparse(bars(:, [1,1,2,2]), ones(size(F))*[1,2,1,2], [Fvec, -Fvec], N, 2));    
    P = P + deltat*Ftot;                            %   Update node positions

    %   Bring outside vertices back to the boundary
    dist  = drectangle(P, L, W);                    %   Find distances
    ind   = dist>0;                                 %   Find vertices outside (d>0)
    dgradx = (drectangle([P(ind, 1)+deps, P(ind, 2)], L, W)-dist(ind))/deps;   %   Numerical
    dgrady = (drectangle([P(ind, 1), P(ind, 2)+deps], L, W)-dist(ind))/deps;   %   gradient
    dgrad2 = dgradx.^2 + dgrady.^2;
    P(ind, :) = P(ind, :) - ...
        [dist(ind).*dgradx./dgrad2, dist(ind).*dgrady./dgrad2];     %   Projection

    %   Termination criterion: All interior nodes move less than dptol (scaled)
    factor = max(sqrt(sum(deltat*Ftot(dist<-geps,:).^2, 2))/h0);    

    if factor<dptol, break; end
    if count>2^13 break; end;
end
%close(h);  

cla; patch('vertices', P, 'faces',t,'edgecol','k','facecol',[.8, .9, 1]); view(0, 90); axis equal;
xlabel('x'); ylabel('y');
count
Triangles = size(t, 1)
quality = min(simpqual(P, t))



function [P, t, I] = plate_refine(P, t, I, error, percentage, par, method)
%   SYNTAX
%   [P, t, I] = plate_refine(P, t, I, error, percentage, par, method)
%   DESCRIPTION
%   This function refines a mesh for a base plate (or other geometry)
%   according to the local error - performs subdivision of triangles with
%   the larger error. Laplacian smoothing of the resulting mesh is used.
%   Inputs:
%   P - initial array of nodes
%   t - initial array of vertices
%   I - initial number of boundary nodes placed first
%   error - array of local errors for every triangle
%   percentage - percentage of triangles to be refined
%   if par = 1 then all edges of the triangle in question are refined
%   if par ~=1 then only edges adjacent to two triangles in question are
%   refined
%   method - method/type of Laplacian smoothing (1, 2, 3, 4)
%   Outputs: refined mesh P, t, and number of boundary nodes I placed first
%   in the list of P
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
    
    [dummy, ind]        = sort(abs(error));   
    ind = ind(round(end*(1-percentage/100))+1:end);                 %   Triangles to be refined
    
    dt     = triangulation(t(:, 1:3), P);                           %   Standard form
    edgesa = edges(dt);                                             %   All edges
    edgesa = sort(edgesa, 2);
    edges0 = freeBoundary(dt);                                      %   Boundary edges    
    edges0 = sort(edges0, 2);    
    edgesi = setxor(edgesa, edges0, 'rows');                        %   Inner edges

    %   Boundary nodes to be added
    Attb   = cell2mat(edgeAttachments(dt, edges0));                 %   Triangles attached to all boundary edges
    edgesr = edges0(ismember(Attb, ind), :);                        %   Indexes into boundary edges to be refined
    nodesb = 0.5*(P(edgesr(:, 1),:) + P(edgesr(:, 2),:));           %   Boundary nodes to be addded 
    nodesb1 = nodesb;
    
    %   Inner nodes to be added
    Attb   = cell2mat(edgeAttachments(dt, edgesi));                 %   Triangles attached to all inner edges
    if par
        edgesr = edgesi(ismember(Attb(:, 1), ind)|...               %   All edges of the triangle in question are refined
                        ismember(Attb(:, 2), ind), :);              %   Indexes into inner edges to be refined
    else
         edgesr = edgesi(ismember(Attb(:, 1), ind)&...              %   Only edges adjacent to two triangles in question are refined
         ismember(Attb(:, 2), ind), :);                             %   Indexes into inner edges to be refined
    end
    nodesb2 = 0.5*(P(edgesr(:, 1),:) + P(edgesr(:, 2),:));          %   Inner nodes to be addded 
    temp = unique(nodesb2, 'rows');
   
    P = [nodesb1; P; nodesb2];
    %   Create the mesh
    DT  = delaunayTriangulation(P);
    t   =  DT.ConnectivityList;
    I = I + size(nodesb, 1);
        
    %   Laplacian smoothing
    M = 20;
    L = max(P(:, 1)) - min(P(:, 1));
    W = max(P(:, 2)) - min(P(:, 2));
  
    geps    = 1e-12*max(L, W);
    quality = zeros(M+1, 1);
    for m = 1:M        
        [P, t]  = laplace_mod(P, I, method); 
        dist    = drectangle(P, L, W);      %  A column of signed distances
        if sum(dist>geps)
            P   = P(dist<geps, :);
            DT  = delaunayTriangulation(P);
            t   =  DT.ConnectivityList;
        end
        quality(m) = min(simpqual(P, t));
        if  m > 1 
            if quality(m) < quality(m-1)
                break;
            end
        end
    end
end
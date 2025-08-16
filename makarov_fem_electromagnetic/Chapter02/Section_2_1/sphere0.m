function [P, t] = sphere0(r, F);
%   SYNTAX
%   P, t] = sphere0(r, F)
%   DESCRIPTION
%   This function generates a sphere mesh from a basic octahedron mesh. [P,
%   t] = sphere0(r, F) creates a sphere mesh of radius R, which has an
%   approximate number of faces F. This function subdivides the basic
%   octahedron so that it has an approximate number of faces F, and shift
%   the vertices to the surface of the boundary of the sphere. Originally
%   written by Mr. Aung Thu Htet.
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
    warning off;
    % Vertices of an octahedron
    P = [ 0,  0 ,  r;
          0,  0 , -r;
          0,  r ,  0;
          0, -r ,  0;
          r,  0 ,  0;
         -r,  0 ,  0; ];
    % Triangulation
    d = DelaunayTri(P);
    [t1, X1] = freeBoundary(d); 
    if (F < 12)
        P = X1;   
        t = t1; 
        Faces = length(t);
        return;
    end
    % Determine the number of loops from faces 
     if (F < 12)
          w = 0;      
     elseif  (F > 12)  && (F <= 48)
          w = 1; 
     elseif  (F > 48)  && (F <= 192)
          w = 2; 
     elseif  (F > 192) && (F <= 768)
          w = 3; 
     elseif  (F > 768) && (F <= 3072)
          w = 4; 
     elseif  (F > 3072) && (F <= 12288)
          w = 5; 
     else w = 6;
     end
    % Loop for dividing edges 
    k = 0; 
    while(k < w)
        % free boundary edges
        E = [];
        for i =1:length(t1);
        E = [E; t1(i,1), t1(i,2) ; t1(i,1), t1(i,3); t1(i,2), t1(i,3) ];
        end
        % free boundary edges  
        % divide edges
        X2 = [];
        for j = 1:length(E);
            X2 = [X2; (X1(E(j,1),1) + X1(E(j,2),1))/2,(X1(E(j,1),2) + X1(E(j,2),2))/2 ,(X1(E(j,1),3) + X1(E(j,2),3))/2 ];
        end
        X2 = [X1;X2];
        % divide edges
        d2 = DelaunayTri(X2);
        [t2, X2] = freeBoundary(d2); 
        t1 = t2;
        X1= X2;  
        k = k + 1;     
    end
    % Convert Cartesian to spherical and shift nodes to sphere boundary
    Xsph = zeros(length(X2),3);
    for i = 1:length(X2);
        xs=X2(i,1);
        ys=X2(i,2);
        zs=X2(i,3);
        azS = atan2(ys,xs);
        elS = atan2(zs,sqrt(xs*xs + ys*ys));
        Xsph(i,1) = r ;
        Xsph(i,2) = azS ;
        Xsph(i,3) = elS ;
    end
    % Convert shifted nodes from spherical to Cartesian
    Xcart = zeros(length(Xsph),3);
    for i = 1:length(Xsph);
        azC= Xsph(i,2);
        elC= Xsph(i,3);
        xc = r .* cos(elC) .* cos(azC);
        yc = r .* cos(elC) .* sin(azC);
        zc = r .* sin(elC);
        Xcart(i,1) = xc;
        Xcart(i,2) = yc;
        Xcart(i,3) = zc;
    end
    % Delaunay and plot final vertices of sphere 
    d3 = DelaunayTri(Xcart);
    [t3, X3] = freeBoundary(d3); 
    % Outputs variables
    P = X3;
    t = t3;
    Faces = length(t3);
end


 
 

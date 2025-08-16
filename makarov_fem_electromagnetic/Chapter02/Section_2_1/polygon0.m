function [P, t] = polygon0(L, W, l, w, indicator)
%   SYNTAX
%   [P, t] = polygon0(L, W, l, w, indicator)
%   DESCRIPTION
%   This function creates a structured centered mesh for a polygon - a blade dipole on the size L by W
%   The feed length is l
%   The feed width is w
%
%   To display the mesh use: fv.faces = t; fv.vertices = P;
%   patch(fv, 'FaceColor', 'y'); axis equal; view(160, 60); grid on;
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
   
    %   Create centered mesh in the xy-plane
    Nx = ceil(2*L/l);   if  mod(Nx,2) Nx = Nx+1; end;
    Ny = ceil(W/w);     if ~mod(Ny,2) Ny = Ny+1; end;
    Nx
    Ny
    if indicator == 0
        x = [0:Nx]/Nx*L;    %   uniform grid
        y = [0:Ny]/Ny*W;    %   uniform grid
    else
        x = L/2*(1 - cos(pi*[0:Nx]/Nx));    %   non-uniform grid
        y = W/2*(1 - cos(pi*[0:Ny]/Ny));    %   non-uniform grid
    end    
    x = x - mean(x);
    y = y - mean(y);
    x(Nx/2)         = -l/2; 
    x(Nx/2+1)       = 0; 
    x(Nx/2+2)       = +l/2;
    y((Ny+1)/2)     = -w/2; 
    y((Ny+1)/2+1)   = +w/2;
    
    P(:, 1) = repmat(x, 1, length(y))';
    P(:, 2) = reshape(repmat(y, length(x), 1), 1, length(x)*length(y))';
    P(:, 1) = P(:, 1) - mean(P(:, 1));
    P(:, 2) = P(:, 2) - mean(P(:, 2));
    P(:, 3) = 0;
    %   Define triangles 
    t = []; temp = [];
    for n = 1:Ny
        for m = 1:Nx
            IND = ((m==Nx/2)|(m==Nx/2+1))&(n~=(Ny+1)/2);
            if IND 
                continue; 
            end            
            if m<Nx/2+1
                if n>Ny/2
                    t1 = [m+1 m         m+Nx+1]' + (n-1)*(Nx+1);
                    t2 = [m+1 m+Nx+2    m+Nx+1]' + (n-1)*(Nx+1);    
                else
                    t1 = [m m+1    m+Nx+2]' + (n-1)*(Nx+1);
                    t2 = [m m+Nx+1 m+Nx+2]' + (n-1)*(Nx+1);   
                end
            else
                if n>Ny/2
                    t1 = [m m+1    m+Nx+2]' + (n-1)*(Nx+1);
                    t2 = [m m+Nx+1 m+Nx+2]' + (n-1)*(Nx+1);
                else
                    t1 = [m+1 m         m+Nx+1]' + (n-1)*(Nx+1);
                    t2 = [m+1 m+Nx+2    m+Nx+1]' + (n-1)*(Nx+1);                
                end
            end                  
            temp = [temp t1 t2];
        end
    end
    t = temp'; 
end

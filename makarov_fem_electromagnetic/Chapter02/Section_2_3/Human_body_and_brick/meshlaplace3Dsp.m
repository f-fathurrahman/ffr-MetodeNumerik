function [Pnew, tnew] = meshlaplace3Dsp(P, t, alpha, beta, iter)
%   SYNTAX 
%   [Pnew, tnew] = meshlaplace3Dsp(P, t, alpha, beta, iter)
%   DESCRIPTION 
%   This function implements a 3D surface-preserving Laplacian smoothing of
%   a triangular mesh based on existing Delaunay connectivity
%   Inputs:
%   P - array of vertices; t - array of faces; alpha, beta, iter - algorithm parameters
%   Outputs:
%   Pnew - new nodal points; tnew = t - array of faces
%   The smoothing method used is given in:
%   "Improved Laplacian Smoothing of Noisy Surface Meshes," J. Vollmer, R. Mencl, and H. Müller
%   EUROGRAPHICS ’99 / P. Brunet and R. Scopigno (Guest Editors) Volume 18 (1999)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    N = size(P, 1);
    %   Triangles attached to every vertex (cell array)
    [si] = meshconnvt(t);
    %   Algorithm
    O = P;
    B = P;
    Iter = iter;
    for iter = 1:Iter
        %iter
        Pprev = P;    
        for m = 1:N                
            index = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
            index(index==m) = [];           %   Indexes into vertices connected to vertex m    
            P(m, :) = sum(Pprev(index, :))/length(index);
            B(m, :) =  P(m, :) - (alpha*O(m, :) + (1-alpha)*Pprev(m, :));        
        end
        for m = 1:N                
            index = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
            index(index==m) = [];           %   Indexes into vertices connected to vertex m               
            P(m, :) = P(m, :) - (beta*B(m, :) + (1-beta)*sum(B(index, :))/length(index));        
        end            
    end
    Pnew = P;
    tnew = t;
end


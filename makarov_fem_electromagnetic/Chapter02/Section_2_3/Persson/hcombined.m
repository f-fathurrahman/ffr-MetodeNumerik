function h = hcombined(P, comb, par)
    %   Compute element size function for a combined structure
    %   0<par<1
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    d = dcombined(P, comb);
    h = 1/par - (1/par-1)*d/mean(abs(d));    %   Relative edge length
end

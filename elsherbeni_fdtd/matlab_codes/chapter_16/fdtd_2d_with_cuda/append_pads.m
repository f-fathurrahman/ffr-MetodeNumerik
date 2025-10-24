function [outarr] = append_pads(inarr, xsize, ysize, xshift, yshift)
outarr = zeros(xsize, ysize);
sz = size(inarr);
outarr(1+xshift:sz(1)+xshift, 1+yshift:sz(2)+yshift) = inarr;
